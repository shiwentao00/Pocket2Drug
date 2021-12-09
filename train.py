import argparse
import yaml
import os
import random
import torch
import torch.nn as nn
from torch.nn.utils.rnn import pad_sequence
from torch.nn.utils.rnn import pack_padded_sequence
from torch.optim.lr_scheduler import ReduceLROnPlateau

from dataloader import pocket_single_loader_gen
from model import Pocket2Drug


def get_args():
    parser = argparse.ArgumentParser("python")

    parser.add_argument("-val_fold",
                        required=False,
                        default=0,
                        help="which fold used for validation")

    return parser.parse_args()


def read_folds(val_fold, data_dir="./data/folds/"):
    # train folds
    folds = list(range(10))
    folds.pop(val_fold)

    # put the data in train folds together
    train_dict = {}
    for fold in folds:
        with open(data_dir + 'pockets_fold{}.yaml'.format(fold), 'r') as f:
            fold_dict = yaml.full_load(f)
        train_dict.update(fold_dict)
    
    # load the validation dict
    with open(data_dir + 'pockets_fold{}.yaml'.format(val_fold), 'r') as f:
        val_dict = yaml.full_load(f)

    return train_dict, val_dict
    

if __name__ == "__main__":
    args = get_args()
    val_fold = int(args.val_fold)
    assert val_fold in list(range(10))
    print('training for cross-validation, validation fold {}.'.format(val_fold))

    # load configuration file
    config_dir = "./train.yaml"
    with open(config_dir, 'r') as f:
        config = yaml.full_load(f)

    # directory for results
    out_dir = config['out_dir']
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    print('results saved in {}.'.format(out_dir))
    trained_model_dir = out_dir + 'trained_model.pt'

    # detect cpu or gpu
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print('device: ', device)

    random.seed(config['seed'])

    # save the configuration file for future reference
    with open(out_dir + 'config.yaml', 'w') as f:
        yaml.dump(config, f)

    # training data files
    pocket_dir = config['pocket_dir']
    pop_dir = config['pop_dir']
    features_to_use = config['features_to_use']

    # load the pocket-smiles pairs
    smiles_train_dict, smiles_val_dict = read_folds(val_fold=val_fold)

    # dataloaders
    batch_size = config['batch_size']
    num_workers = os.cpu_count()
    num_workers = int(min(batch_size, num_workers))
    print('number of workers to load data: ', num_workers)
    trainloader, train_size = pocket_single_loader_gen(
        smiles_train_dict,
        pocket_dir,
        pop_dir,
        features_to_use,
        vocab=config['vocab'],
        vocab_path=config['vocab_path'],
        batch_size=batch_size, shuffle=False,
        num_workers=num_workers
    )
    print('size of train set: ', train_size)

    valloader, val_size = pocket_single_loader_gen(
        smiles_val_dict,
        pocket_dir,
        pop_dir,
        features_to_use,
        vocab=config['vocab'],
        vocab_path=config['vocab_path'],
        batch_size=batch_size, shuffle=False,
        num_workers=num_workers
    )
    print('size of val set: ', val_size)

    # model initialization
    encoder_config = config['encoder_config']
    decoder_config = config['decoder_config']
    model = Pocket2Drug(encoder_config, decoder_config).to(device)

    # load pretrained encoder
    if encoder_config['pretrain']:
        print('loading pretrained GNN as encoder...')
        loaded_gnn = torch.load(
            encoder_config['pretrained_model'],
            map_location=torch.device(device)
        )
        load_info = model.load_state_dict(loaded_gnn, strict=False)
        print(load_info)
        print('Pretrained GNN for encoder is loaded.')
    else:
        print('No pretraining for encoder GNN.')

    # load pretrained decoder
    if decoder_config['pretrain']:
        print('loading pretrained RNN as decoder...')
        model.decoder.load_state_dict(
            torch.load(
                decoder_config['pretrained_model'],
                map_location=torch.device(device)
            )
        )
        print('Pretrained RNN for decoder is loaded.')
    else:
        print('No pretraining for decoder RNN')

    # the optimizer
    learning_rate = config['learning_rate']
    weight_decay = config['weight_decay']
    optimizer = torch.optim.Adam(
        model.parameters(),
        lr=learning_rate,
        weight_decay=weight_decay,
        amsgrad=True
    )

    # the learning rate scheduler
    scheduler = ReduceLROnPlateau(
        optimizer,
        mode='min',
        factor=0.5,
        patience=3,
        cooldown=10,
        min_lr=0.0001,
        verbose=True
    )

    loss_function = nn.CrossEntropyLoss(reduction='sum')

    # get the index of padding
    PADDING_IDX = config['decoder_config']['num_embeddings'] - 1

    # train and validation, the results are saved.
    train_losses = []
    val_losses = []
    best_val_loss, best_val_epoch = float('inf'), None
    num_epoch = config['num_epoch']
    print('begin training...')
    for epoch in range(1, 1 + num_epoch):
        # train
        model.train()
        train_loss = 0
        for data in trainloader:
            optimizer.zero_grad()
            data = data.to(device)
            smiles = data.y

            # the lengths are decreased by 1 because we don't
            # use <eos> for input and we don't need <sos> for
            # output during traning.
            lengths = [len(x) - 1 for x in smiles]

            # pad the sequences
            smiles = [torch.tensor(x) for x in smiles]
            smiles = pad_sequence(
                smiles, batch_first=True,
                padding_value=PADDING_IDX
            ).to(device)

            # forward
            preds = model(data, smiles, lengths)

            # The <sos> token is removed before packing, because
            # we don't need <sos> of output during training.
            # Note that the lengths are already decreased by 1.
            targets = pack_padded_sequence(
                smiles[:, 1:],
                lengths,
                batch_first=True,
                enforce_sorted=False
            ).data

            loss = loss_function(preds, targets)
            loss.backward()
            optimizer.step()
            train_loss += loss.item()  # * data.num_graphs
        train_losses.append(train_loss / train_size)

        # validation
        model.eval()
        val_loss = 0
        for data in valloader:
            data = data.to(device)
            smiles = data.y

            # the lengths are decreased by 1 because we don't
            # use <eos> for input and we don't need <sos> for
            # output during traning.
            lengths = [len(x) - 1 for x in smiles]

            # pad the sequences
            smiles = [torch.tensor(x) for x in smiles]
            smiles = pad_sequence(
                smiles, batch_first=True,
                padding_value=PADDING_IDX
            ).to(device)

            # forward
            preds = model(data, smiles, lengths)

            # The <sos> token is removed before packing, because
            # we don't need <sos> of output during training.
            # Note that the lengths are already decreased by 1.
            targets = pack_padded_sequence(
                smiles[:, 1:],
                lengths,
                batch_first=True,
                enforce_sorted=False
            ).data

            loss = loss_function(preds, targets)
            val_loss += loss.item()  # * data.num_graphs
        val_losses.append(val_loss / val_size)

        print('epoch {}, train loss: {}, val loss: {}.'.format(
            epoch, train_losses[-1], val_losses[-1]))

        # update the saved model upon best validation loss
        if val_losses[-1] <= best_val_loss:
            best_val_epoch = epoch
            best_val_loss = val_losses[-1]
            torch.save(model.state_dict(), trained_model_dir)
            print('model saved at epoch {}'.format(epoch))

        scheduler.step(val_losses[-1])

    # save train and validation losses
    loss_history = [train_losses, val_losses]
    with open(out_dir + 'loss.yaml', 'w') as f:
        yaml.dump(loss_history, f)

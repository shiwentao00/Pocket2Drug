import yaml
import os
import random
import torch
import torch.nn as nn
from torch.nn.utils.rnn import pad_sequence
from torch.nn.utils.rnn import pack_padded_sequence
from torch.optim.lr_scheduler import ReduceLROnPlateau
from tqdm import tqdm

from dataloader import pocket_loader_gen
from model import Pocket2Drug


if __name__ == "__main__":
    # detect cpu or gpu
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print('device: ', device)

    config_dir = "./train.yaml"
    with open(config_dir, 'r') as f:
        config = yaml.full_load(f)

    random.seed(config['seed'])

    # directory for results
    out_dir = config['out_dir']
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    trained_model_dir = out_dir + 'trained_model.pt'

    # save the configuration file for future reference
    with open(out_dir + 'config.yaml', 'w') as f:
        yaml.dump(config, f)

    # training data files
    pocket_dir = config['pocket_dir']
    pop_dir = config['pop_dir']
    features_to_use = config['features_to_use']

    # load the pocket-smiles pairs
    smiles_dir = config['smiles_dir']
    with open(smiles_dir, 'r') as f:
        smiles_dict = yaml.full_load(f)

    # exclude pockets used in case study
    excluded_pockets = config['excluded_pockets']
    with open(excluded_pockets, 'r') as f:
        excluded_pockets = yaml.full_load(f)
    for pocket in excluded_pockets:
        smiles_dict.pop(pocket)

    # dataloaders
    batch_size = config['batch_size']
    num_workers = os.cpu_count()
    num_workers = int(min(batch_size, num_workers))
    print('number of workers to load data: ', num_workers)
    trainloader, valloader, train_size, val_size = pocket_loader_gen(
        smiles_dict,
        pocket_dir,
        pop_dir,
        features_to_use,
        vocab=config['vocab'],
        vocab_path=config['vocab_path'],
        batch_size=batch_size, shuffle=False,
        test_split=config['test_split'],
        num_workers=num_workers
    )

    # model and training configuration
    encoder_config = config['encoder_config']
    decoder_config = config['decoder_config']
    model = Pocket2Drug(encoder_config, decoder_config).to(device)
    learning_rate = config['learning_rate']
    weight_decay = config['weight_decay']
    loss_function = nn.CrossEntropyLoss(reduction='sum')

    # the optimizer
    optimizer = torch.optim.Adam(
        model.parameters(),
        lr=learning_rate,
        weight_decay=weight_decay,
        amsgrad=True
    )
    # for name, param in model.named_parameters():
    #    if param.requires_grad:
    #        # print(name, param.data)
    #        print(name, param)

    # the learning rate scheduler
    scheduler = ReduceLROnPlateau(
        optimizer,
        mode='min',
        factor=0.5,
        patience=5,
        cooldown=20,
        min_lr=0.0001,
        verbose=True
    )

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
        for data in tqdm(trainloader):
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

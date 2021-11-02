import sys
import os
import time
import argparse
import numpy as np
import tensorflow as tf
import keras
from keras.optimizers import Adam
from keras.layers import Conv1D
from keras.layers import GlobalAveragePooling1D
from keras.models import Model
from keras.layers import Dense, Input, Dropout
from keras.layers.normalization import BatchNormalization
from keras import optimizers
from keras import regularizers
from keras.regularizers import l2
from keras.callbacks import TensorBoard
from keras.callbacks import EarlyStopping, ModelCheckpoint, ReduceLROnPlateau


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dataset', type=str,
	help="npz file with input data")
parser.add_argument('-g','--gpu', default=-1, type=int,
	help='the GPU number, -1 indicates CPU')
parser.add_argument('-p','--prefix', default='folder_train', type=str,
	help='prefix of the output folder where the models are saved')

args = parser.parse_args()

gpu = args.gpu

if gpu > -1:
	device = '/gpu:%i' % gpu
	os.environ['CUDA_VISIBLE_DEVICES']=str(args.gpu)
else:
	device = '/cpu:0'

#os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

if args.prefix == None:
	parser.print_help()
	sys.stderr.write("Please specify the dataset npz file!\n")
	sys.exit(1)


timestr = time.strftime("%Y%m%d-%H%M%S")
logdir = '%s_%s/' % (args.prefix, timestr)
os.mkdir(logdir)


def build_model():
    input_seq = Input(shape=(600, 4))
    features = Conv1D(128, 10, padding='same', kernel_initializer="he_normal",
                   kernel_regularizer=l2(1e-4), activation='relu')(input_seq)
    features = BatchNormalization()(features)
    features = Dropout(0.1)(features)
    layer = GlobalAveragePooling1D()(features)
    layer = Dense(128, activation='relu')(layer)
    layer = BatchNormalization()(layer)
    layer = Dropout(0.1)(layer)
    layer = Dense(2, activation='relu')(layer)
    layer = BatchNormalization()(layer)
    layer = Dropout(0.1)(layer)
    layer = Dense(1, activation='sigmoid')(layer)
    model = Model(inputs=input_seq, outputs=layer)
    model.summary()
    adam = Adam(lr=0.0001)
    model.compile(optimizer=adam, loss='binary_crossentropy', metrics=['accuracy'])
    return model


data = np.load(args.dataset, allow_pickle=True, encoding='bytes')
Y = data['cl']
folds = data['partition']
partitions_interval = np.arange(6)
inner_partitions_interval = partitions_interval[partitions_interval != 5]
# Inner cross-validation
for val_partition in inner_partitions_interval: 
	os.mkdir('%spartition_%i' % (logdir, val_partition))
	model_file = "%spartition_%i/model" % (logdir, val_partition)
	model_file2 = "%spartition_%i/model_nomgpu" % (logdir, val_partition)
	log_file = "%spartition_%i/" % (logdir, val_partition)
	train_partition = inner_partitions_interval[inner_partitions_interval != val_partition]
	train_set = np.in1d(folds.ravel(), train_partition).reshape(folds.shape)
	val_set = np.where(folds == val_partition)

	# Load training data
	X_train = data['X'][train_set]
	y_train = Y[train_set]

	# Load validation data
	X_val = data['X'][val_set]
	y_val = Y[val_set]
    
	model = build_model()
	checkpointer = ModelCheckpoint(filepath=model_file,verbose=1, save_best_only=True)

	earlyStopping = keras.callbacks.EarlyStopping(monitor='val_loss',
                                              patience=15, verbose=1,
                                              mode='auto')
	lr_reducer = keras.callbacks.ReduceLROnPlateau(monitor='val_loss',
                                               factor=0.1, patience=5,
                                               verbose=0, mode='auto',
                                               min_delta=0.0001,
                                               cooldown=0, min_lr=0)
	tensorboard = TensorBoard(log_dir=log_file,
                          histogram_freq=0,
                          write_graph=True,
                          write_images=False)
	history = model.fit(X_train, y_train, validation_data=(X_val, y_val), epochs=150, verbose=1, batch_size=64,
						callbacks=[checkpointer, earlyStopping, tensorboard, lr_reducer]).history

model.save(model_file2)
print("Saved model to disk")
    

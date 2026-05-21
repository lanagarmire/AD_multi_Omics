import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from keras.layers import Input, Dense
from keras.models import Model
from keras.optimizers import Adam
from keras.metrics import MeanSquaredError
from keras import backend as K

from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from scipy import stats

bottleneck_dim = 5000


# ---- load data ----
data = pd.read_csv('/nfs/turbo/umms-lgarmire/home/yhdu/Bowei_NAS/az_omics_new/az_omics/mrna_data.csv')

data_t = data.T
data_t.columns = data_t.iloc[0, :]
data_t = data_t.iloc[1:, :]
data_mat = data_t.to_numpy().astype('float32')

input_dim = data_mat.shape[1]

x_train, x_test = train_test_split(data_mat, test_size=0.2, random_state=42)
print("Train/Test:", x_train.shape, x_test.shape)

def clear_sess():
    K.clear_session()
    import gc
    gc.collect()

clear_sess()

# ---- model ----
input_data = Input(shape=(input_dim,))

#encoded = Dense(1024, activation='linear', name='encoder1')(input_data)
#enc = Dense(512, activation='linear', name='encoder2')(encoded)
#enc1 = Dense(bottleneck_dim, activation='linear', name='encoder3')(enc)
#dec1 = Dense(512, activation='linear', name='decoder3')(enc1)
#dec = Dense(1024, activation='linear', name='decoder2')(dec1)
#decoded = Dense(input_dim, activation='linear', name='decoder1')(dec)

# Simple model: 20032 → 1024 → 20032
enc1 = Dense(bottleneck_dim, activation='linear', name='encoder_simple')(input_data)
decoded = Dense(input_dim, activation='linear', name='decoder_simple')(enc1)

autoencoder = Model(inputs=input_data, outputs=decoded)

autoencoder.compile(
    optimizer=Adam(learning_rate=0.001),
    loss='mse',
    metrics=[MeanSquaredError()]
)

autoencoder.summary()

history = autoencoder.fit(
    x_train, x_train,
    epochs=10,
    batch_size=32,
    validation_data=(x_test, x_test),
    verbose=1
)

# ---- save history ----
loss = history.history['loss']
val_loss = history.history['val_loss']
epochs = np.arange(1, len(loss) + 1)

np.savez(
    f"autoencoder_history_{bottleneck_dim}.npz",
    loss=np.array(loss),
    val_loss=np.array(val_loss),
    epochs=epochs
)

# ---- plotting ----
plt.figure(figsize=(7, 5))
plt.plot(epochs, loss, label='Training Loss', marker='o')
plt.plot(epochs, val_loss, label='Validation Loss', marker='o')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.title(f'Autoencoder Training Convergence, Bottleneck = {bottleneck_dim}')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(
    f"/nfs/turbo/umms-lgarmire/home/yhdu/Bowei_NAS/EFIGA/Review/Autoencoder/convergence_{bottleneck_dim}.png",
    dpi=300
)
plt.close()

# ---- reconstruct & R2/PCC ----
x_decoded = autoencoder.predict(x_test)

ov_r2 = []
n_eval = min(len(x_test), len(x_decoded))
for i in range(n_eval):
    x_original = x_test[i]
    x_recons = x_decoded[i]
    r2 = r2_score(x_original.tolist(), x_recons.tolist())
    ov_r2.append(r2)

print("overall R2", sum(ov_r2) / n_eval)

x_original = x_test[0]
x_recons = x_decoded[0]
pcc, _ = stats.pearsonr(x_original.tolist(), x_recons.tolist())
print("PCC first sample:", pcc)

# ---- reduced dimensions (metagenes) ----
layer_name_1 = 'encoder_simple'
intermediate_layer_model_1 = Model(
    inputs=autoencoder.input,
    outputs=autoencoder.get_layer(layer_name_1).output
)

intermediate_output_1 = intermediate_layer_model_1(data_mat.astype('float32'))
intermediate_output_1 = intermediate_output_1.numpy()

metagene = [f'metagene_{k}' for k in range(1, bottleneck_dim + 1)]

data_hmrna = pd.DataFrame(intermediate_output_1, columns=metagene)
data_hmrna = data_hmrna.loc[:, (data_hmrna != 0).any(axis=0)]

data_hmrna.to_csv(f'Autoencorder_metagene_681p_{bottleneck_dim}hf.csv', index=False)

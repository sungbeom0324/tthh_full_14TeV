# AE without label (unsupervised), reconstructing MNEST Images of numbers 0~9.
# AE is learning topological common features of 0~9.
import numpy as np
import tensorflow as tf
from tensorflow.keras.layers import Input, Dense
from tensorflow.keras.models import Model
from tensorflow.keras.datasets import mnist
import matplotlib.pyplot as plt

(x_train, _), (x_test, _) = mnist.load_data()

# Data normalization
x_train = x_train.astype('float32') / 255.
x_test = x_test.astype('float32') / 255.

# Convert image into vector
x_train = x_train.reshape((len(x_train), np.prod(x_train.shape[1:])))
x_test = x_test.reshape((len(x_test), np.prod(x_test.shape[1:])))

# Input shape
input_img = Input(shape=(784,))
# Latent space
encoded = Dense(3, activation='relu')(input_img)
# Reconstrtion
decoded = Dense(784, activation='sigmoid')(encoded)

# MODEL : AE, EN, DE #
# AE(I/O)
autoencoder = Model(input_img, decoded)

# EN
encoder = Model(input_img, encoded)

# place holder for the input of decoder
encoded_input = Input(shape=(3,))
# place holder for the output of decoder
decoder_layer = autoencoder.layers[-1]
# DE
decoder = Model(encoded_input, decoder_layer(encoded_input))

# Compile Model Set
autoencoder.compile(optimizer='adam', loss='binary_crossentropy')
autoencoder.fit(x_train, x_train,
                epochs=50,
                batch_size=256,
                shuffle=True,
                validation_data=(x_test, x_test))

# Save model
autoencoder.save('my_model.h5')

# Prediction #
encoded_imgs = encoder.predict(x_test)
decoded_imgs = decoder.predict(encoded_imgs)

# RE #
reconstruction_error = np.mean(np.square(x_test - decoded_imgs))
print("Reconstruction Error (MSE):", reconstruction_error)

n = 10  # of numbers to display
plt.figure(figsize=(20, 4))
for i in range(n):
    # Original Image
    ax = plt.subplot(2, n, i + 1)
    plt.imshow(x_test[i].reshape(28, 28))
    plt.gray()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    # Reconstructed Image
    ax = plt.subplot(2, n, i + n + 1)
    plt.imshow(decoded_imgs[i].reshape(28, 28))
    plt.gray()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
plt.savefig("fig0106.pdf")


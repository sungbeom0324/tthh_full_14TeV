# AE reconstructing MNIST images of numbers 0~9.  
# It is unsupervised though, extract dixcriminative regions.
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow.keras.layers import Input, Dense
from tensorflow.keras.models import Model

# 데이터 생성
(x_train, y_train), (x_test, y_test) = tf.keras.datasets.mnist.load_data()
x_train, x_test = x_train / 255.0, x_test / 255.0  # 이미지 스케일링

# Autoencoder 모델 정의
input_img = Input(shape=(784,))
encoded = Dense(128, activation='relu')(input_img)
encoded = Dense(64, activation='relu')(encoded)
encoded = Dense(3, activation='relu')(encoded)  # 3차원 latent space

decoded = Dense(64, activation='relu')(encoded)
decoded = Dense(128, activation='relu')(decoded)
decoded = Dense(784, activation='sigmoid')(decoded)

autoencoder = Model(input_img, decoded)

# 모델 컴파일
autoencoder.compile(optimizer='adam', loss='binary_crossentropy')

# 모델 학습
autoencoder.fit(x_train.reshape((len(x_train), np.prod(x_train.shape[1:]))),
                x_train.reshape((len(x_train), np.prod(x_train.shape[1:]))),
                epochs=50,
                batch_size=256,
                shuffle=True,
                validation_data=(x_test.reshape((len(x_test), np.prod(x_test.shape[1:]))),
                                 x_test.reshape((len(x_test), np.prod(x_test.shape[1:])))))

# Encoder 모델 정의 (3차원 latent space 출력)
encoder = Model(input_img, encoded)

# 잠재 공간에서 각 숫자를 나타내는 색상 부여
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple', 'brown']

# 테스트 데이터를 잠재 공간에 매핑하여 3D 시각화
encoded_imgs = encoder.predict(x_test.reshape((len(x_test), np.prod(x_test.shape[1:]))))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(10):
    indices = np.where(y_test == i)
    ax.scatter(encoded_imgs[indices, 0], encoded_imgs[indices, 1], encoded_imgs[indices, 2], c=colors[i], label=str(i))

ax.legend()
plt.savefig("fig_0107.pdf")


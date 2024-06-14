import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, MinMaxScaler
from sklearn.metrics import confusion_matrix
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.utils import to_categorical
import matplotlib.pyplot as plt

# 读取和预处理数据
cts_gbm = pd.read_csv("bioinfoHW/DataFiles/TCGA-GBM.htseq_counts.tsv", sep="\t", index_col="Ensembl_ID")
cts_gbm = np.round(2**cts_gbm - 1)

cts_pcpg = pd.read_csv("bioinfoHW/DataFiles/TCGA-PCPG.htseq_counts.tsv", sep="\t", index_col="Ensembl_ID")
cts_pcpg = np.round(2**cts_pcpg - 1)

# 随机选取表达量不为0的基因作为识别标准
np.random.seed(123)
positive_rows = cts_gbm[(cts_gbm != 0).all(axis=1) & (cts_pcpg != 0).all(axis=1)].index
sampled_rows = np.random.choice(positive_rows, size=1000, replace=False)

# 建立标签列表
def label_data(cts, label):
    thecolnames = cts.columns
    thelabel = []
    for var in thecolnames:
        if var[13] == "0":
            thelabel.append("tumor")
        else:
            thelabel.append("normal")
    coldata = pd.DataFrame(thelabel, index=thecolnames, columns=["condition"])
    coldata["condition"] = coldata["condition"].apply(lambda x: label if x == "tumor" else 0)
    return coldata

coldata_gbm = label_data(cts_gbm, 1)
coldata_pcpg = label_data(cts_pcpg, 2)

# 筛选数据和归一化
def normalize_data(cts, coldata, sampled_rows):
    cts_sampled = cts.loc[sampled_rows]
    scaler = MinMaxScaler()
    cts_sampled_norm = scaler.fit_transform(cts_sampled.T).T
    cts_sampled_norm_df = pd.DataFrame(cts_sampled_norm, index=sampled_rows, columns=cts.columns)
    new_coldata = coldata.join(cts_sampled_norm_df.T)
    idx_train, idx_valid = train_test_split(new_coldata.index, train_size=0.8, random_state=123)
    return new_coldata.loc[idx_train], new_coldata.loc[idx_valid]

data_gbm_train, data_gbm_valid = normalize_data(cts_gbm, coldata_gbm, sampled_rows)
data_pcpg_train, data_pcpg_valid = normalize_data(cts_pcpg, coldata_pcpg, sampled_rows)

# 合并训练和验证数据
cts_combined_train = pd.concat([data_gbm_train, data_pcpg_train])
cts_combined_valid = pd.concat([data_gbm_valid, data_pcpg_valid])
le = LabelEncoder()
cts_combined_train['condition'] = le.fit_transform(cts_combined_train['condition'])
cts_combined_valid['condition'] = le.transform(cts_combined_valid['condition'])

# 准备数据
train_x = cts_combined_train.drop(columns='condition').values
train_y = to_categorical(cts_combined_train['condition'], num_classes=3)
valid_x = cts_combined_valid.drop(columns='condition').values
valid_y = to_categorical(cts_combined_valid['condition'], num_classes=3)

# 构建全连接神经网络模型
model = Sequential([
    Dense(512, activation='relu', input_shape=(train_x.shape[1],)),
    Dropout(0.5),
    Dense(256, activation='relu'),
    Dropout(0.5),
    Dense(128, activation='relu'),
    Dropout(0.5),
    Dense(3, activation='softmax')
])

# 编译模型
model.compile(loss='categorical_crossentropy',
              optimizer='adam',
              metrics=['accuracy'])

# 训练模型
history = model.fit(train_x, train_y,
                    epochs=30,
                    batch_size=32,
                    validation_data=(valid_x, valid_y))

# 绘制训练历史
plt.plot(history.history['accuracy'], label='accuracy')
plt.plot(history.history['val_accuracy'], label='val_accuracy')
plt.xlabel('Epoch')
plt.ylabel('Accuracy')
plt.ylim([0, 1])
plt.legend(loc='lower right')
plt.show()

# 模型评估
test_loss, test_acc = model.evaluate(valid_x, valid_y, verbose=2)
print('\nTest accuracy:', test_acc)

# 输出混淆矩阵
predictions = np.argmax(model.predict(valid_x), axis=1)
cm = confusion_matrix(cts_combined_valid['condition'], predictions)
print(cm)

# 保存模型
model.save("fcnn_model.h5")

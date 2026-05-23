"""
Лабораторная работа №6: Распознавание жестов (LeapGestRecog)
Классификация 10 жестов с использованием MobileNetV2 и аугментации.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow.keras.applications import MobileNetV2
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, GlobalAveragePooling2D
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint, ReduceLROnPlateau
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import train_test_split
import seaborn as sns
from PIL import Image
import glob

# ==================== ПАРАМЕТРЫ ====================
DATA_DIR = "./leapGestRecog"          # путь к датасету
IMG_SIZE = (224, 224)                 # размер для MobileNetV2
BATCH_SIZE = 32
EPOCHS = 30
VALIDATION_SPLIT = 0.2
USE_PRETRAINED = True
RANDOM_SEED = 42

# ==================== ЗАГРУЗКА ДАННЫХ ====================
print("🔍 Сбор путей к изображениям...")
all_images = glob.glob(os.path.join(DATA_DIR, "*", "*", "*.png"))
print(f"Найдено всего изображений: {len(all_images)}")

def extract_class_name(path):
    """Извлекает имя жеста из пути: .../01_palm/... → palm"""
    dir_name = os.path.basename(os.path.dirname(path))  # "01_palm"
    if "_" in dir_name:
        parts = dir_name.split("_", 1)
        if len(parts) == 2:
            return parts[1]   # "palm"
    return dir_name

class_names = sorted(set(extract_class_name(p) for p in all_images))
num_classes = len(class_names)
print("Обнаруженные классы:", class_names)
print(f"Всего классов: {num_classes}")

# Создаём список (путь, метка)
data = [(img, extract_class_name(img)) for img in all_images]
label_to_id = {name: i for i, name in enumerate(class_names)}

X_paths = [item[0] for item in data]
y = [label_to_id[item[1]] for item in data]

# Разделение на train/val
X_train_paths, X_val_paths, y_train, y_val = train_test_split(
    X_paths, y, test_size=VALIDATION_SPLIT, stratify=y, random_state=RANDOM_SEED
)
print(f"Train: {len(X_train_paths)}, Val: {len(X_val_paths)}")

# Функция загрузки и предобработки
def load_and_preprocess(path, label, augment=False):
    image = tf.io.read_file(path)
    image = tf.image.decode_png(image, channels=3)
    image = tf.image.resize(image, IMG_SIZE)
    image = tf.cast(image, tf.float32) / 255.0
    if augment:
        image = tf.image.random_flip_left_right(image)
        image = tf.image.random_brightness(image, 0.2)
        image = tf.image.random_contrast(image, 0.8, 1.2)
        # Простой поворот на -15..15 градусов
        angle = tf.random.uniform([], -0.2, 0.2)
        image = tf.image.rot90(image, k=tf.cast(angle * 10, tf.int32) % 4)
    return image, label

def make_dataset(paths, labels, batch_size, augment=False):
    dataset = tf.data.Dataset.from_tensor_slices((paths, labels))
    dataset = dataset.map(lambda p, l: load_and_preprocess(p, l, augment),
                          num_parallel_calls=tf.data.AUTOTUNE)
    dataset = dataset.batch(batch_size).prefetch(tf.data.AUTOTUNE)
    return dataset

train_ds = make_dataset(X_train_paths, y_train, BATCH_SIZE, augment=True)
val_ds = make_dataset(X_val_paths, y_val, BATCH_SIZE, augment=False)

# ==================== МОДЕЛЬ ====================
if USE_PRETRAINED:
    print("\n🤖 Используем предобученную MobileNetV2...")
    base = MobileNetV2(input_shape=(*IMG_SIZE, 3), include_top=False, weights='imagenet')
    base.trainable = False
    model = Sequential([
        base,
        GlobalAveragePooling2D(),
        Dense(128, activation='relu'),
        Dropout(0.5),
        Dense(num_classes, activation='softmax')
    ])
else:
    print("\n🏗️ Своя CNN...")
    model = Sequential([
        tf.keras.layers.Conv2D(32, (3,3), activation='relu', input_shape=(*IMG_SIZE, 3)),
        tf.keras.layers.MaxPooling2D(2,2),
        tf.keras.layers.Conv2D(64, (3,3), activation='relu'),
        tf.keras.layers.MaxPooling2D(2,2),
        tf.keras.layers.Conv2D(128, (3,3), activation='relu'),
        tf.keras.layers.MaxPooling2D(2,2),
        tf.keras.layers.Flatten(),
        tf.keras.layers.Dense(256, activation='relu'),
        tf.keras.layers.Dropout(0.5),
        Dense(num_classes, activation='softmax')
    ])

model.compile(optimizer=tf.keras.optimizers.Adam(0.001),
              loss='sparse_categorical_crossentropy',
              metrics=['accuracy'])
model.summary()

# ==================== КОЛБЭКИ ====================
callbacks = [
    EarlyStopping(monitor='val_accuracy', patience=5, restore_best_weights=True, verbose=1),
    ModelCheckpoint('best_gesture_model.keras', monitor='val_accuracy', save_best_only=True, verbose=1),
    ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=3, verbose=1)
]

# ==================== ОБУЧЕНИЕ ====================
print("\n🚀 Начинаем обучение...")
history = model.fit(train_ds, validation_data=val_ds, epochs=EPOCHS, callbacks=callbacks, verbose=1)

# Загрузка лучшей модели
model.load_weights('best_gesture_model.keras')

# ==================== ОЦЕНКА ====================
print("\n📊 Оценка на валидации...")
val_loss, val_acc = model.evaluate(val_ds, verbose=0)
print(f"Validation accuracy: {val_acc*100:.2f}%")

# Предсказания для отчёта
y_true, y_pred = [], []
for imgs, lbls in val_ds:
    preds = model.predict(imgs, verbose=0)
    y_true.extend(lbls.numpy())
    y_pred.extend(np.argmax(preds, axis=1))

print("\n📈 Classification Report:")
print(classification_report(y_true, y_pred, target_names=class_names))

# Матрица ошибок
cm = confusion_matrix(y_true, y_pred)
plt.figure(figsize=(10,8))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=class_names, yticklabels=class_names)
plt.title('Confusion Matrix')
plt.ylabel('True')
plt.xlabel('Predicted')
plt.savefig('confusion_matrix.png', dpi=150)
plt.close()

# Графики обучения
plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
plt.plot(history.history['accuracy'], label='Train Acc')
plt.plot(history.history['val_accuracy'], label='Val Acc')
plt.legend(); plt.grid(); plt.title('Accuracy')
plt.subplot(1,2,2)
plt.plot(history.history['loss'], label='Train Loss')
plt.plot(history.history['val_loss'], label='Val Loss')
plt.legend(); plt.grid(); plt.title('Loss')
plt.tight_layout()
plt.savefig('training_history.png', dpi=150)
plt.close()
print("✅ Графики сохранены: training_history.png, confusion_matrix.png")

# ==================== ПРЕДСКАЗАНИЕ СВОИХ ФОТО ====================
def predict_gesture(image_path, model, class_names, img_size=IMG_SIZE):
    if not os.path.exists(image_path):
        print(f"❌ {image_path} не найден")
        return
    try:
        img = Image.open(image_path).convert('RGB')
        img = img.resize(img_size)
        arr = np.array(img) / 255.0
        arr = np.expand_dims(arr, axis=0)
        probs = model.predict(arr, verbose=0)[0]
        idx = np.argmax(probs)
        print(f"🖐️ {os.path.basename(image_path)} → {class_names[idx]} (уверенность: {probs[idx]*100:.1f}%)")
        # Топ-3
        top3 = np.argsort(probs)[-3:][::-1]
        print(f"   Топ-3: {', '.join([f'{class_names[i]} ({probs[i]*100:.1f}%)' for i in top3])}")
    except Exception as e:
        print(f"Ошибка при {image_path}: {e}")

print("\n" + "="*60)
print("ПРЕДСКАЗАНИЕ ДЛЯ ВАШИХ ФОТОГРАФИЙ (my_hand_*.jpg)")
print("="*60)

gesture_files = sorted(glob.glob("my_hand_*.jpg"))
if gesture_files:
    for f in gesture_files:
        predict_gesture(f, model, class_names)
else:
    print("Файлы my_hand_*.jpg не найдены.")
    print("Положите свои фото в текущую папку и назовите my_hand_1.jpg, my_hand_2.jpg и т.д.")

print("\n🏁 Готово! Модель сохранена как best_gesture_model.keras")

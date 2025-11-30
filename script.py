import pandas as pd      # бібліотека для аналізу та маніпуляцій даних
import numpy as np       # для генерації чисел і масивів
import matplotlib.pyplot as plt  # для побудови графіків
import ROOT              # для створення ROOT-файлів і TTree
from array import array  # для передачі масивів у ROOT TTree

# -------------------------------
# Параметри генерації
# -------------------------------
np.random.seed(42)
means = [1, 2, 5]  # очікувані значення для експоненційного розподілу
n_samples = 1000   # кількість значень у кожному розподілі

# -------------------------------
# Генерація даних та створення DataFrame (Pandas)
# -------------------------------
data_dict = {f'mean_{mean}': np.random.exponential(scale=mean, size=n_samples)
             for mean in means}   # створюємо словник зі списками значень
df = pd.DataFrame(data_dict)       # перетворюємо словник у DataFrame

# -------------------------------
# Базова інспекція даних (Pandas)
# -------------------------------
print("Перші 5 рядків DataFrame:")
print(df.head())           # перегляд перших рядків DataFrame
print("\nОписова статистика:")
print(df.describe())       # описова статистика: mean, std, min, max, percentiles

# -------------------------------
# Додавання нового стовпця — трансформація даних (Pandas)
# -------------------------------
df['sum_means'] = df.sum(axis=1)  # сума всіх значень у рядку
print("\nПерші 5 рядків з новим стовпцем 'sum_means':")
print(df.head())

# -------------------------------
# Групування та бінування даних (Pandas)
# -------------------------------
df['sum_bin'] = pd.cut(df['sum_means'], bins=5)  # бінування значень суми
grouped = df.groupby('sum_bin').mean()           # середнє по бінованих групах
print("\nСереднє по групах бінів 'sum_means':")
print(grouped)

# -------------------------------
# Збереження DataFrame у CSV (Pandas)
# -------------------------------
df.to_csv("exponential_data.csv", index=False)
print("\nDataFrame збережено у 'exponential_data.csv'.")

# -------------------------------
# Побудова гістограм (matplotlib)
# -------------------------------
plt.figure(figsize=(10,6))

for mean in means:
    plt.hist(df[f'mean_{mean}'], bins=30, alpha=0.5, label=f'Mean={mean}')

plt.title("Експоненційні розподіли для різних значень очікуваної величини")
plt.xlabel("Значення")
plt.ylabel("Частота")
plt.legend()
plt.grid(True)

plt.savefig("exponential_histograms.png", dpi=300)  # збереження графіку
plt.show()
print("Гістограма збережена у 'exponential_histograms.png'.")

# -------------------------------
# Збереження даних у ROOT TTree (ROOT)
# -------------------------------
root_file = ROOT.TFile("exponential_distributions.root", "RECREATE")  # створення ROOT-файлу
tree = ROOT.TTree("tree", "Exponential Distributions")                 # створення TTree

arrays = {mean: array('f', df[f'mean_{mean}']) for mean in means}     # перетворюємо Pandas Series у масиви для ROOT
for mean in means:
    tree.Branch(f'mean_{mean}', arrays[mean], f'mean_{mean}[{n_samples}]/F')

tree.Fill()        # заповнюємо TTree даними
tree.Write()       # запис у ROOT-файл
root_file.Close()  # закриваємо файл
print("ROOT-файл 'exponential_distributions.root' створено успішно.")


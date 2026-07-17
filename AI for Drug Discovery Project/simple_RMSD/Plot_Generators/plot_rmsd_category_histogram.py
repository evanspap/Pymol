"""Plot RMSD category distribution from CSV data."""

import csv
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt


def load_categories(csv_file):
    categories = []
    with open(csv_file, newline='') as csvfile:
        # Use a permissive encoding fallback for files with non-UTF-8 bytes
        content = csvfile.read().replace('\r\n', '\n')
    for line in content.split('\n'):
        if not line or line.startswith('#'):
            continue
        if line.startswith('TargetFrame,Frame,RMSD,Category'):
            continue
        parts = line.split(',')
        if len(parts) < 4:
            continue
        categories.append(parts[3].strip())
    return categories


def plot_category_distribution(categories, output_path):
    counts = Counter(categories)
    category_order = ['<=2 Å', '2-3 Å', '3-4 Å', '4-5 Å', '>5 Å']
    category_labels = [label for label in category_order if label in counts]
    values = [counts[label] for label in category_labels]

    plt.figure(figsize=(10, 6))
    bars = plt.bar(category_labels, values, color='tab:blue', edgecolor='black')

    for bar in bars:
        height = bar.get_height()
        plt.annotate(f'{height:,}',
                     xy=(bar.get_x() + bar.get_width() / 2, height),
                     xytext=(0, 6),
                     textcoords='offset points',
                     ha='center', va='bottom', fontsize=10)

    plt.title('RMSD Category Distribution')
    plt.xlabel('RMSD Category')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def main():
    base_dir = Path(__file__).parent
    csv_file = base_dir / 'rmsd_all_vs_all.csv'
    if not csv_file.exists():
        raise FileNotFoundError(f'CSV file not found: {csv_file}')

    categories = load_categories(csv_file)
    if not categories:
        raise ValueError('No categories found in CSV file.')

    output_image = base_dir / 'rmsd_category_histogram.png'
    plot_category_distribution(categories, output_image)
    print(f'Saved histogram to: {output_image}')


if __name__ == '__main__':
    main()

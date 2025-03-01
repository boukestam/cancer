import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import seaborn as sns
import os
import umap
from sklearn.metrics import silhouette_score
import scipy.cluster.hierarchy as sch
from sklearn.cluster import AgglomerativeClustering

data_dir = "data"
result_dir = "results"

# Load the activities file
df = pd.read_csv(os.path.join(data_dir, 'all_25_activities.txt'), sep='\t', index_col=0)

# Calculate the total mutations per sample and normalize
df['Total_Mutations'] = df.sum(axis=1)
df = df.div(df['Total_Mutations'], axis=0)
df = df.drop(columns=['Total_Mutations'])

# Load the metadata file
meta_df = pd.read_csv(os.path.join(result_dir, 'metadata.csv'))
df = df.merge(meta_df, left_on='Samples', right_on='Tumor_Sample_Barcode')

print(df.head())

# Plotting
df_grouped = df.drop(columns=['Tumor_Sample_Barcode']).groupby('Project_Code').mean()
df_grouped.plot(kind='bar', stacked=True, figsize=(12, 8))
plt.title('Average Signature Contributions per Cancer Type')
plt.xlabel('Cancer Type')
plt.ylabel('Proportion of Mutations')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(result_dir, "signature_contributions.png"), dpi=300, bbox_inches="tight")
plt.clf()

# Scale data
numeric_df = df.drop(columns=['Tumor_Sample_Barcode', 'Project_Code'])
scaler = StandardScaler()
df_scaled = scaler.fit_transform(numeric_df)

# Perform K-Means clustering
# kmeans = KMeans(n_clusters=20, random_state=42)  # Adjust the number of clusters as needed
# numeric_df['Cluster'] = kmeans.fit_predict(df_scaled)
# df['Cluster'] = kmeans.labels_

# Elbow method

# inertia = []
# K_range = range(2, 20)

# for k in K_range:
#   kmeans_test = KMeans(n_clusters=k, random_state=42, n_init=10)
#   kmeans_test.fit(df_scaled)
#   inertia.append(kmeans_test.inertia_)

# plt.plot(K_range, inertia, marker='o')
# plt.xlabel('Number of Clusters')
# plt.ylabel('Inertia')
# plt.title('Elbow Method for Optimal K')
# plt.savefig(os.path.join(result_dir, "elbow.png"), dpi=300, bbox_inches="tight")
# plt.clf()

# Hierarchical clustering
linkage_matrix = sch.linkage(df_scaled, method='ward')
hierarchical = AgglomerativeClustering(n_clusters=10, metric='euclidean', linkage='ward')
numeric_df['Cluster'] = hierarchical.fit_predict(df_scaled)
df['Cluster'] = hierarchical.labels_

# Silhouette score
score = silhouette_score(df_scaled, hierarchical.labels_)
print(f"Silhouette Score: {score:.2f}")

# Visualize the clustering with a heatmap
plt.figure(figsize=(12, 8))
sns.heatmap(numeric_df.sort_values('Cluster'), cmap='viridis', cbar=True)
plt.title('Heatmap of Signature Contributions with Clustering')
plt.xlabel('Signatures')
plt.ylabel('Samples')
plt.savefig(os.path.join(result_dir, "cluster_heatmap.png"), dpi=300, bbox_inches="tight")
plt.clf()

# Count number of samples per (Cluster, Cancer Type)
cluster_counts = df.groupby(['Cluster', 'Project_Code']).size().reset_index(name='Count')

# Pivot to make it more readable
cluster_pivot = cluster_counts.pivot(index='Cluster', columns='Project_Code', values='Count').fillna(0)

# Display as a heatmap
plt.figure(figsize=(12, 6))
sns.heatmap(cluster_pivot, annot=True, fmt=".0f", cmap="Blues", linewidths=0.5)
plt.title("Cancer Type Distribution Across Clusters")
plt.xlabel("Cancer Type (Project_Code)")
plt.ylabel("Cluster")
plt.xticks(rotation=45, ha='right')
plt.savefig(os.path.join(result_dir, "cluster_cancer_distribution.png"), dpi=300, bbox_inches="tight")
plt.clf()

plt.figure(figsize=(12, 6))
sns.barplot(data=cluster_counts, x='Cluster', y='Count', hue='Project_Code')
plt.title("Cancer Type Distribution in Clusters")
plt.xlabel("Cluster")
plt.ylabel("Number of Samples")
plt.legend(title="Cancer Type", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig(os.path.join(result_dir, "cluster_cancer_barplot.png"), dpi=300, bbox_inches="tight")
plt.clf()

# Reduce to 2D using UMAP
umap_model = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
df_umap = umap_model.fit_transform(df_scaled)

# Plot clusters
plt.figure(figsize=(10, 6))
scatter = plt.scatter(df_umap[:,0], df_umap[:,1], c=df['Cluster'], cmap='tab10', alpha=0.6)
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.title("UMAP Visualization of Mutational Signatures Clusters")
plt.colorbar(scatter, label="Cluster")
plt.savefig(os.path.join(result_dir, "cluster_umap.png"), dpi=300, bbox_inches="tight")
plt.clf()
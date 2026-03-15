## Optional features

### Network hard limit bypass

Backend : On the build_knowledge_graph.py 

```python
# Change this:
# de_genes = get_top_de_genes(study_id, limit=60)

# To this:
de_genes = get_top_de_genes(study_id, limit=2000) # Actually build a large network
```
Frontend : 6_Knowledge_Graph.py

```python
# Add this under your Edge sliders
max_nodes = st.slider("Max Visible Genes (Nodes)", min_value=20, max_value=500, value=60)

# Modify the node processing loop
nodes = []
for n in nodes_data[:max_nodes]: # <--- Slice the nodes here dynamically
    logfc = n.get("logFC", 0)
    color = "#ff4b4b" if logfc > 0 else "#2b83ff"
    
    nodes.append(Node(
        id=n["id"],
        label=n["label"],
        size=25,
        color=color,
        title=f"Gene: {n['label']}\nLogFC: {logfc}"
    ))
```
This will allow you build a larger network 

### Rank Based / Quantile Normalization

This script is designed to take the massive MIDBASE master matrix and a new user's sample, and mathematically force the user's sample into the same technical space without altering the biological signature.



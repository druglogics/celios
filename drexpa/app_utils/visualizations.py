"""
Visualization utilities for DREXPA results
"""

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def plot_targets_per_drug(drug_targets_df):
    """
    Create a bar chart showing number of targets per drug.

    Args:
        drug_targets_df: DataFrame with drug_name and targets columns

    Returns:
        plotly.graph_objects.Figure
    """
    if drug_targets_df.empty:
        return go.Figure().add_annotation(text="No data available")

    # Count targets per drug
    if "drug_name" in drug_targets_df.columns and "targets" in drug_targets_df.columns:
        # Parse targets (assuming comma-separated or list)
        target_counts = []

        for _, row in drug_targets_df.iterrows():
            drug = row["drug_name"]
            targets = row["targets"]

            # Handle different formats
            if isinstance(targets, str):
                target_list = [t.strip() for t in targets.split(",") if t.strip()]
            elif isinstance(targets, list):
                target_list = targets
            else:
                target_list = []

            target_counts.append({"drug": drug, "n_targets": len(target_list)})

        counts_df = pd.DataFrame(target_counts).sort_values("n_targets", ascending=False)

        fig = px.bar(
            counts_df,
            x="drug",
            y="n_targets",
            title="Number of Targets per Drug",
            labels={"drug": "Drug", "n_targets": "Target Count"},
            color="n_targets",
            color_continuous_scale="Viridis",
        )
        fig.update_layout(
            xaxis_tickangle=-45,
            height=500,
            hovermode="x unified"
        )
        return fig

    return go.Figure().add_annotation(text="Required columns not found")


def plot_drug_target_heatmap(drug_targets_df):
    """
    Create a heatmap of drug-target interactions.

    Args:
        drug_targets_df: DataFrame with drug_name and targets columns

    Returns:
        plotly.graph_objects.Figure
    """
    if drug_targets_df.empty:
        return go.Figure().add_annotation(text="No data available")

    if "drug_name" not in drug_targets_df.columns or "targets" not in drug_targets_df.columns:
        return go.Figure().add_annotation(text="Required columns not found")

    # Build drug-target matrix
    drug_targets = {}
    all_targets = set()

    for _, row in drug_targets_df.iterrows():
        drug = row["drug_name"]
        targets = row["targets"]

        # Parse targets
        if isinstance(targets, str):
            target_list = [t.strip() for t in targets.split(",") if t.strip()]
        elif isinstance(targets, list):
            target_list = targets
        else:
            target_list = []

        drug_targets[drug] = set(target_list)
        all_targets.update(target_list)

    # Limit to top targets for readability
    all_targets = sorted(list(all_targets))[:30]  # Top 30 targets

    # Create matrix
    matrix = []
    drugs = list(drug_targets.keys())

    for target in all_targets:
        row = [1 if target in drug_targets.get(drug, set()) else 0 for drug in drugs]
        matrix.append(row)

    # Create heatmap
    fig = go.Figure(
        data=go.Heatmap(
            z=matrix,
            x=drugs,
            y=all_targets,
            colorscale="RdYlGn",
            colorbar=dict(title="Target<br>Present"),
        )
    )

    fig.update_layout(
        title="Drug-Target Interaction Matrix",
        xaxis_title="Drug",
        yaxis_title="Target",
        height=600,
        xaxis_tickangle=-45,
    )

    return fig


def plot_network_diagram(drug_targets_df):
    """
    Create a network diagram of drug-target interactions.

    Args:
        drug_targets_df: DataFrame with drug_name and targets columns

    Returns:
        plotly.graph_objects.Figure
    """
    if drug_targets_df.empty:
        return go.Figure().add_annotation(text="No data available")

    if "drug_name" not in drug_targets_df.columns or "targets" not in drug_targets_df.columns:
        return go.Figure().add_annotation(text="Required columns not found")

    # Build network
    nodes = []
    edges_x = []
    edges_y = []
    node_colors = []
    node_sizes = []
    edge_trace_data = []

    drug_positions = {}
    target_positions = {}

    # Add drug nodes
    for idx, drug in enumerate(drug_targets_df["drug_name"].unique()):
        nodes.append(f"Drug: {drug}")
        drug_positions[drug] = idx
        node_colors.append("blue")
        node_sizes.append(15)

    # Add target nodes and edges
    target_idx = len(drug_positions)
    all_targets = set()

    for _, row in drug_targets_df.iterrows():
        drug = row["drug_name"]
        targets = row["targets"]

        # Parse targets
        if isinstance(targets, str):
            target_list = [t.strip() for t in targets.split(",") if t.strip()]
        elif isinstance(targets, list):
            target_list = targets
        else:
            target_list = []

        for target in target_list[:5]:  # Limit to top 5 targets per drug for clarity
            all_targets.add(target)

            # Add target node if not already added
            if target not in target_positions:
                target_positions[target] = target_idx
                nodes.append(f"Target: {target}")
                node_colors.append("red")
                node_sizes.append(10)
                target_idx += 1

            # Add edge
            edge_trace_data.append({
                "drug": drug,
                "target": target,
                "drug_pos": drug_positions[drug],
                "target_pos": target_positions[target]
            })

    # Create a simple layout using plotly scatter
    # Position drugs on the left, targets on the right
    node_x = []
    node_y = []
    node_text = []

    import math

    # Position drugs
    n_drugs = len(drug_positions)
    for drug, idx in drug_positions.items():
        angle = (idx / max(n_drugs, 1)) * 2 * math.pi
        x = 1 + 0.3 * math.cos(angle)
        y = 0.3 * math.sin(angle)
        node_x.append(x)
        node_y.append(y)
        node_text.append(drug)

    # Position targets
    n_targets = len(target_positions)
    for target, idx in target_positions.items():
        angle = (idx / max(n_targets, 1)) * 2 * math.pi
        x = -1 + 0.3 * math.cos(angle)
        y = 0.3 * math.sin(angle)
        node_x.append(x)
        node_y.append(y)
        node_text.append(target)

    # Create edges
    edge_x = []
    edge_y = []
    for edge in edge_trace_data:
        x0 = node_x[edge["drug_pos"]]
        y0 = node_y[edge["drug_pos"]]
        x1 = node_x[n_drugs + edge["target_pos"]]
        y1 = node_y[n_drugs + edge["target_pos"]]

        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    # Create figure
    fig = go.Figure()

    # Add edges
    fig.add_trace(
        go.Scatter(
            x=edge_x,
            y=edge_y,
            mode="lines",
            line=dict(width=0.5, color="gray"),
            hoverinfo="none",
            showlegend=False,
        )
    )

    # Add nodes
    fig.add_trace(
        go.Scatter(
            x=node_x,
            y=node_y,
            mode="markers+text",
            text=node_text,
            textposition="top center",
            hoverinfo="text",
            marker=dict(
                size=node_sizes,
                color=node_colors,
                opacity=0.8,
            ),
            showlegend=False,
        )
    )

    fig.update_layout(
        title="Drug-Target Interaction Network",
        showlegend=False,
        hovermode="closest",
        margin=dict(b=20, l=5, r=5, t=40),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        height=600,
    )

    return fig

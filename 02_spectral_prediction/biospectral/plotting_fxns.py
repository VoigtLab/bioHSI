import matplotlib.pyplot as plt
import plotly.graph_objects as go


def plot_spectra(spectra_df, save_name=None, **kwargs):
    fig = plt.figure()
    for cpd in spectra_df.columns:
        plt.plot(x, spectra_df[cpd], kwargs)
    plt.ylabel(r"$\epsilon$", fontsize=16)
    plt.xlabel("Wavelength (nm)", fontsize=16)
    plt.tight_layout()
    plt.legend()

    if save_name:
        plt.savefig(save_name, dpi=600)
    plt.show()


def plot_spectra_interactive(spectra_df, save_name=None, **kwargs):
    fig = go.Figure()

    # Add traces
    for cpd in spectra_df.columns:
        fig.add_trace(
            go.Scatter(
                kwargs,
                x=spectra_df.index,
                y=spectra_df[cpd],
                mode="lines",
                name=cpd,
                text=cpd,
            )
        )
    fig.update_layout(showlegend=False)
    return fig

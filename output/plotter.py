"""

"""
import matplotlib.pyplot as plt
import plotly.graph_objects as go


class Plotter:

    @staticmethod
    def plot_pulse_input_and_output(xy_par,
                                    graph_title: str, x_graph: str,
                                    y_graph: str):
        """

        Args:
            xy_par:
            y_graph:
            x_graph:
            graph_title:

        Returns:

        """
        fig = go.Figure()

        for i in xy_par:
            fig.add_trace(go.Scatter(
                x=i[0],
                y=i[1],
                name=i[2],
            ))

        fig.update_layout(title=graph_title,
                          xaxis_title=x_graph,
                          yaxis_title=y_graph)

        fig.show()
        fig.write_html(graph_title + ".html")

    @staticmethod
    def plot_power_output_and_delta_output(z, power_output, delta_output, name):
        """

        Args:
            z:
            power_output:
            delta_output:
            name:

        Returns:

        """
        fig, ax1 = plt.subplots()
        color = 'red'
        ax1.set_xlabel('distance (km)')
        ax1.set_ylabel('Output Power', color=color)
        ax1.plot(z, power_output, color=color)
        ax1.tick_params(axis='y', labelcolor=color)

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        color = 'blue'
        ax2.set_ylabel('Delta', color=color)  # we already handled the x-label with ax1
        ax2.plot(z, delta_output, color=color)
        ax2.tick_params(axis='y', labelcolor=color)

        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        fig.savefig(name)

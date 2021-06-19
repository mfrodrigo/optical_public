"""

"""
import matplotlib.pyplot as plt


class Plotter:

    @staticmethod
    def plot_pulse_input_and_output(t, pulse_input, pulse_output, name):
        """

        Args:
            t:
            pulse_input:
            pulse_output:
            name:

        Returns:

        """
        fig = plt.figure()
        plt.plot(t, abs(pulse_input[:, 0]) ** 2, label='Input')
        plt.plot(t, abs(pulse_output[:, 0]) ** 2, label='Output')
        plt.title('Gaussian Pulse ')
        plt.xlabel(r'$(t-\beta_1z)/T_0$')
        plt.ylabel('|u1(z,t)|^2/P_0')
        plt.legend()
        plt.grid(True)
        fig.savefig(name)

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

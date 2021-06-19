"""

"""
import pandas as pd


class Tables:

    @staticmethod
    def export_table_results(distance, power_output, delta_output, name):
        data = pd.DataFrame({'Distância': distance,
                             'Potência de saída': power_output,
                             'Largura a meia altura': delta_output})
        data.to_csv(name, decimal=",")

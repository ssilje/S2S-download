# environment dependencies
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

# local dependencies
from S2S.local_configuration_H import config
from .create_domain_file import make_dir

def save_fig(filename):

    make_dir(config['SAVEFIG'])
    plt.savefig(config['SAVEFIG']+filename)
    plt.close()
    print('Figure stored at: '+config['SAVEFIG']+filename)

def qq_plot(df,figure_name):

    sns.set_theme(style="ticks")
    sns.lmplot(x="observations", y="model", data=df,
                height=2, aspect=1.5, col='season',row='lead_time',
                scatter_kws={"s": 50, "alpha": 0.3})
    save_fig(figure_name+'_qq.png')

def line_plot(df,figure_name,var_name='values'):

    sns.set_theme(style="ticks")
    sns.relplot(
        data=df,
        x="lead_time", y=var_name, col="month",
        kind="line", palette="crest", linewidth=4, zorder=5,
        col_wrap=3, height=2, aspect=1.5)
    save_fig(figure_name+'_line.png')

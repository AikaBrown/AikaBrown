{
  "nbformat": 4,
  "nbformat_minor": 0,
  "metadata": {
    "colab": {
      "name": "bioreactor dynamics",
      "provenance": [],
      "collapsed_sections": [],
      "authorship_tag": "ABX9TyO+YZpDFmhXuSsjIYkWbTf5",
      "include_colab_link": true
    },
    "kernelspec": {
      "name": "python3",
      "display_name": "Python 3"
    },
    "language_info": {
      "name": "python"
    }
  },
  "cells": [
    {
      "cell_type": "markdown",
      "metadata": {
        "id": "view-in-github",
        "colab_type": "text"
      },
      "source": [
        "<a href=\"https://colab.research.google.com/github/AikaBrown/AikaBrown/blob/main/bioreactor_dynamics.py\" target=\"_parent\"><img src=\"https://colab.research.google.com/assets/colab-badge.svg\" alt=\"Open In Colab\"/></a>"
      ]
    },
    {
      "cell_type": "code",
      "execution_count": null,
      "metadata": {
        "id": "9lTaBNE_JIH8"
      },
      "outputs": [],
      "source": [
        "def RungeKutta4(T, X, S, y, h, μmax, Ks, V, F, Sr, I):\n",
        "  \"\"\"\n",
        "  Description:\n",
        "    In order to know the dynamic behaviour of a biochemical reactor, we have to simulate the developed model. \n",
        "    The model structure of the CONTINUOUS STIRRED TANK BIOREACTOR (CSTB) is comprised of the following set of equations:\n",
        "\n",
        "    1) Biomass balance -->    dX/dT = x*(μ-D)\n",
        "    2) Substrate balance -->  dS/dT = D*(Sf-s)- μ*x/y\n",
        "    μ = μmax*s/(km+s)\n",
        "\n",
        "    For solving the above differential–algebraic system, we have to used the fourth-order Runge–Kutta method. \n",
        "    This fuction takes the arguments presented in Args section and solve the fourth-order Runge–Kutta method to obtain a Data Frame\n",
        "    with values of Time (T) in hours, Biomass (X) in g/L and substrate (S) in g/L. \n",
        "\n",
        "\n",
        "  Args:\n",
        "    T:    Time in hours\n",
        "    X:    Initial biomass concentration in g/L\n",
        "    S:    Initial substrate concentration in g/L\n",
        "    y:    The ratio of the amount of biomass produced to the amount of substrate consumed(g biomass/g substrate)\n",
        "    h:    Step size (ΔT) in hours\n",
        "    μmax: The maximum growth rate of the microorganism (h^-1)\n",
        "    Ks:   substrate concentration at which 1/2 μmax is reached (g/L)\n",
        "    V:    Volume of the reactor in L.\n",
        "    F:    Feed flow in L/h\n",
        "    Sr:   Feed substrate concentration (g/L)\n",
        "    I:    Number of approximations\n",
        "\n",
        "  Returns:\n",
        "    DataFrame with values of Time (T) in hours, Biomass (X) in g/L and substrate (S) in g/L\n",
        "\n",
        "  Extra:\n",
        "    this function uses pandas abbreviated as pd:\n",
        "     before use the function, run \"import pandas as pd\"\n",
        "  \"\"\"\n",
        "\n",
        "  D=F/V # Dilution (h^-1)\n",
        "  Valores_prueba = []\n",
        "  Valores= [T, X, S]\n",
        "  Valores_prueba.append(Valores)\n",
        "\n",
        "  for i in range(I):\n",
        "    \n",
        "    #Calculate K1 and L1\n",
        "    f1 = (μmax*X*S)/(Ks+S)-(X*D)\n",
        "    f2 = D*(Sr-S) - (μmax*S/(Ks+S))*(X/y)\n",
        "    K1=h*f1\n",
        "    L1=h*f2\n",
        "\n",
        "    #Calculate K2 and L2\n",
        "    nT = T+h/2\n",
        "    nX = X+K1/2\n",
        "    nS = S+L1/2\n",
        "    nf1 = (μmax*nX*nS)/(Ks+nS)-(nX*D)\n",
        "    nf2 = D*(Sr-nS) - (μmax*nS/(Ks+nS))*(nX/y)\n",
        "    K2=h*nf1\n",
        "    L2=h*nf2\n",
        "    \n",
        "    #Calculate K3 and L3\n",
        "    nT = T+h/2\n",
        "    nX = X+K2/2\n",
        "    nS = S+L2/2\n",
        "    nf1 = (μmax*nX*nS)/(Ks+nS)-(nX*D)\n",
        "    nf2 = D*(Sr-nS) - (μmax*nS/(Ks+nS))*(nX/y)\n",
        "    K3=h*nf1\n",
        "    L3=h*nf2\n",
        "\n",
        "    #Calculate K4 and L4\n",
        "    nT = T+h\n",
        "    nX = X+K3\n",
        "    nS = S+L3\n",
        "    nf1 = (μmax*nX*nS)/(Ks+nS)-(nX*D)\n",
        "    nf2 = D*(Sr-nS) - (μmax*nS/(Ks+nS))*(nX/y)\n",
        "    K4=h*nf1\n",
        "    L4=h*nf2\n",
        "    #Calcular los nuevos valores de T, X, S\n",
        "   \n",
        "    T=T+h\n",
        "    X= X+(1/6)*(K1+2*K2+2*K3+K4)\n",
        "    S= S+(1/6)*(L1+2*L2+2*L3+L4)\n",
        "    Valores= [T, X, S]\n",
        "    Valores_prueba.append(Valores)\n",
        "    df = pd.DataFrame(Valores_prueba, columns=['T', 'X', 'S'])\n",
        "  return(df)    "
      ]
    },
    {
      "cell_type": "code",
      "source": [
        "def biomass_substrate_plot(df):\n",
        "  '''\n",
        "  Description:\n",
        "    Plot Biomass vs substrate concentration\n",
        "  Args:\n",
        "    df: a DataFrame with 3 columns and numerical index: \n",
        "        1) T: Time in hours\n",
        "        2) X: Biomass in g/L\n",
        "        3) S: Substrate in g/L\n",
        "  Returns:\n",
        "    a plot x= Biomass g/L, y= Substrate g/L\n",
        "  \n",
        "  Extra:\n",
        "    this function uses matplotlib.pyplot abbreviated as plt:\n",
        "      before use the function, run \"import matplotlib.pyplot as plt\"\n",
        "  '''\n",
        "  \n",
        "  plt.style.use(\"ggplot\")\n",
        "  fig, ax=plt.subplots()\n",
        "  ax.plot(df['X'], df['S'], linestyle=\"--\", color='b')\n",
        "  ax.set_xlabel(\"Biomasa (g/L)\")\n",
        "  ax.set_ylabel(\"Sustrato (g/L)\")\n",
        "  ax.set_title(\"Biomasa vs Sustrato\\nDiagrama de Fases\")"
      ],
      "metadata": {
        "id": "TIby0BO7JTo0"
      },
      "execution_count": null,
      "outputs": []
    },
    {
      "cell_type": "code",
      "source": [
        "def biomass_time_plot(df):\n",
        "  '''\n",
        "  Description:\n",
        "    Plot Biomass concentration vs time \n",
        "  Args:\n",
        "    df: a DataFrame with 3 columns and numerical index: \n",
        "        1) T: Time in hours\n",
        "        2) X: Biomass in g/L\n",
        "        3) S: Substrate in g/L\n",
        "  Returns:\n",
        "    a plot x= Biomass g/L, y= Substrate g/L\n",
        "  \n",
        "  Extra:\n",
        "    this function uses matplotlib.pyplot abbreviated as plt:\n",
        "      before use the function, run \"import matplotlib.pyplot as plt\"\n",
        "  '''\n",
        "# Graficar biomasa vs tiempo\n",
        "  plt.style.use(\"ggplot\")\n",
        "  fig, ax=plt.subplots()\n",
        "  ax.plot(df['T'], df['S'], linestyle=\"--\", color='y')\n",
        "  ax2 = ax.twinx()\n",
        "  ax2.plot(df['T'], df['X'], linestyle=\"--\", color='g')\n",
        "  ax.set_xlabel(\"Tiempo (h)\")\n",
        "  ax.set_ylabel(\"Sustrato/Biomasa (g/L)\")\n",
        "  ax.set_title(\"Biomasa y Sustrato vs Tiempo \")\n",
        "  plt.show()"
      ],
      "metadata": {
        "id": "be36bFpPJXmI"
      },
      "execution_count": null,
      "outputs": []
    },
    {
      "cell_type": "code",
      "source": [
        ""
      ],
      "metadata": {
        "id": "3PFUG_L_Jadu"
      },
      "execution_count": null,
      "outputs": []
    }
  ]
}
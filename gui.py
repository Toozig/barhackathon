from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import BD as bd
import mcmc_discrete as mcmc

iterations, kT, dT = 0, 0, 0
all_proteins = ["CD28", "CD80", "CD86", "CTLA-4", "PD-1", "PD-L1", "PD-L2", "IL-2RA", "IL-12R",
                "IL-2", "IL-12"]

isValid = False
MCMC = 1
BD = 2


def show_proteins(window, check_buttons, choices_list):
    """
    show protein choices in the given window
    :param window: gui window to update
    :param check_buttons: empty check buttons list to update
    :param choices_list: list of boolean variables for the choices of the check buttons
    """
    proteins_label = Label(window, text="Choose your proteins, and click 'Done':",
                           font=('Calibri bold', 18), background="lavender")
    proteins_label.place(x=300, y=280)
    for i in range(len(all_proteins)):
        c = Checkbutton(window, text=all_proteins[i], var=choices_list[i],
                        font=('Calibri', 14),
                        background="lavender")
        check_buttons.append(c)

    x = 100
    y = 350
    for i in range(len(check_buttons)):
        if x > 900:
            x = 100
            y = y + 40
        check_buttons[i].place(x=x, y=y)
        x += 150


def show_int_args(window):
    iter_label = Label(window, text="Enter number of iterations: ", font=('Calibri bold', 16),
                       background="lavender")
    iter_label.place(x=100, y=100)
    iter_num = Entry(window, bd=2, font=('Calibri', 14), width=6)
    iter_num.place(x=400, y=100)

    kT_label = Label(window, text="Enter kT value: ", font=('Calibri bold', 16),
                     background="lavender")
    kT_label.place(x=100, y=160)
    k = Entry(window, bd=2, font=('Calibri', 14), width=6)
    k.place(x=300, y=160)

    dT_label = Label(window, text="Enter Δt value: ", font=('Calibri bold', 16),
                     background="lavender")
    dT_label.place(x=100, y=220)
    d = Entry(window, bd=2, font=('Calibri', 14), width=6)
    d.place(x=300, y=220)
    return iter_num, k, d


def init_choices(choices_list):
    """
    initialize the given list for the user's choices (through gui window)
    """
    for i in range(len(all_proteins)):
        x = BooleanVar()
        x.set(False)
        choices_list.append(x)


def init_window(window, title, text, x, y):
    """
    initialize the gui window, with the given parameters
    :param window: the gui window
    :param title: title to set
    :param text: text to write at the top of the window (header)
    :param x: the position (row) to place the text(title)
    :param y: the position (column) to place the text(title)
    """
    def on_closing():
        if messagebox.askokcancel("Quit", "Do you want to quit?"):
            window.destroy()
            exit(0)

    window.protocol("WM_DELETE_WINDOW", on_closing)

    window.configure(background="lavender")
    window.geometry("1000x600")
    window.title(title)
    label = Label(window, text=text, font=("Calibri bold", 30),
                  background="lavender")
    label.place(x=x, y=y)


def BDgui(chosen_proteins):
    """
    Graphic user interface for running BD algorithm, with the user's inputs
    :param chosen_proteins: an empty list to be filled with the protein the user chose
    """
    global isValid
    window = Tk()
    choices_list = []
    init_choices(choices_list)
    init_window(window, "BD algorithm", "Choose your BD attributes", 280, 0)
    iter_num, k, d = show_int_args(window)
    check_buttons = []
    show_proteins(window, check_buttons, choices_list)

    def is_done():
        global iterations, kT, dT, isValid
        try:
            iterations = int(iter_num.get())
            kT = float(k.get())
            dT = float(d.get())
            for i in range(len(choices_list)):
                if choices_list[i].get():
                    chosen_proteins.append(all_proteins[i])
            window.destroy()
            isValid = True
        except:
            error_label = Label(window, text="Invalid input, try again", font=('Calibri bold', 16),
                                background="lavender", fg="red")
            error_label.place(x=540, y=530)
            isValid = False

    btn = Button(window, text="Done", bg="black", fg="LavenderBlush2", font=('Calibri bold', 18),
                 command=is_done, background="LavenderBlush4")
    btn.place(x=450, y=520)
    window.mainloop()


def MAINgui():
    """
    Graphic user interface for running BD or MCMC algorithm, according to the user's choice
    """
    window = Tk()
    choose_var = IntVar()
    init_window(window, "Main window", "Choose Algorithm", 330, 0)
    label = Label(window, text="Choose a program to run, and click 'Done':", font=('Calibri bold', 18),
                  background="lavender")
    mcmc = Radiobutton(window, text='MCMC algorithm', value=MCMC, font=('Calibri', 16),
                       background="lavender", variable=choose_var)
    bd = Radiobutton(window, text='BD algorithm', value=BD, font=('Calibri', 16),
                     background="lavender", variable=choose_var)
    label.place(x=280, y=100)
    mcmc.place(x=400, y=180)
    bd.place(x=400, y=240)
    mcmc.deselect()
    bd.deselect()

    def is_done():
        window.destroy()

    btn = Button(window, text="Done", bg="black", fg="LavenderBlush2", font=('Calibri bold', 18),
                 command=is_done, background="LavenderBlush4")
    btn.place(x=450, y=520)
    window.mainloop()
    return choose_var.get()


def MCMCgui(chosen_proteins):
    """
    Graphic user interface for running MCMC algorithm, with the user's inputs
    :param chosen_proteins: an empty list to be filled with the protein the user chose
    """
    global dT
    window = Tk()
    init_window(window, "MCMC algorithm", "Choose your MCMC attributes", 250, 0)
    label = Label(window, text="Choose a csv to upload :", font=('Calibri bold', 18),
                  background="lavender")
    window.filename = ""

    def load_file():
        window.filename = filedialog.askopenfilename(initialdir="/", title="Select file",
                                                     filetypes=(("csv files", "*.csv"), ("all files", "*.*")))
        label = Label(window, text=str(window.filename), font=('Calibri', 13),
                      background="lavender", fg="gray25")
        label.place(x=300, y=150)

    file_but = Button(window, text="load file", font=('Calibri', 14),
                      background="thistle3", command=load_file)
    label.place(x=300, y=100)
    file_but.place(x=610, y=98)
    label = Label(window, text="Enter Δt value: ", font=('Calibri bold', 18),
                  background="lavender")
    d = Entry(window, bd=2, font=('Calibri', 14), width=6)
    label.place(x=300, y=200)
    d.place(x=500, y=205)
    choices_list = []
    init_choices(choices_list)
    check_buttons = []
    show_proteins(window, check_buttons, choices_list)

    def is_done():
        global dT, isValid
        try:
            dt_error, file_err = "Invalid Δt value", "Invalid file: must be csv file"
            dT = float(d.get())
            suffix = ".csv"
            if window.filename == "" or window.filename[-4:] != suffix:
                raise Exception(file_err)
            isValid = True
            for i in range(len(choices_list)):
                if choices_list[i].get():
                    chosen_proteins.append(all_proteins[i])
            window.destroy()
        except Exception as exp:
            isValid = False
            if str(exp) != file_err:
                message = dt_error
            else:
                message = exp
            error_label = Label(window, text=message, font=('Calibri bold', 16),
                                background="lavender", fg="red")
            error_label.place(x=540, y=530)
            isValid = False

    btn = Button(window, text="Done", bg="black", fg="LavenderBlush2", font=('Calibri bold', 18),
                 command=is_done, background="LavenderBlush4")
    btn.place(x=450, y=520)
    window.mainloop()
    return window.filename


def histogram(x):

    x = pd.Series(x)
    x.plot.hist(grid=True, bins=20, rwidth=0.9, color='#607c8e')
    plt.title('Simulated Fret results')
    # plt.xlabel('Δt')
    plt.xlabel('time')
    plt.ylabel('luminous intensity')
    plt.xticks(np.arange(0, iterations * dT, dT))
    plt.grid(axis='y', alpha=0.75)
    plt.show()


def read_excel(protein_vec):
    df = pd.read_csv("newMatrix.csv", index_col=0)
    filtered = df.loc[protein_vec, protein_vec]
    # filtered = filtered.iloc[col_ind]
    filtered['empty state'] = 1
    protein_dict = {}
    for i in range(filtered.shape[0] + 1):
        if i == filtered.shape[0]:
            protein_dict[i] = 'empty state'
        else:
            protein_dict[i] = filtered.columns[i]
    return filtered, protein_dict


if __name__ == '__main__':
    """
    The main function that runs the program.
    At first the main gui window is opened, letting the user to choose the program he wants
    to run (BD or MCMC). Then, according to the user's choice, the matching window is opened
    and the user can enter his arguments.
    Finally, an histogram plot is shown, according to the chosen program's output
    """
    algo = MAINgui()
    chosen_proteins = []
    file_path = ""
    if algo == MCMC:
        run = MCMCgui
    else:
        run = BDgui
    while not isValid:
        chosen_proteins = []
        file_path = run(chosen_proteins)
    if algo == MCMC:
        vec = pd.read_csv(file_path)
        mcmc_obj = mcmc.MCMC(len(chosen_proteins), kT, iterations, chosen_proteins,dT, list(vec))
        results = mcmc_obj.mcmc()
        best_config, hist_arr = results[1], results[3]
        print("Best configurations: ", best_config)
    else:
        matrix, idxs = read_excel(chosen_proteins)
        bd_obj = bd.BD(iterations, kT, dT, (0, 0), len(chosen_proteins), idxs, chosen_proteins, matrix)
        hist_arr = bd_obj.BD_algorithm()[1]
    histogram(hist_arr)





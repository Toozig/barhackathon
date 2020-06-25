from tkinter import *
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import BD

iterations, kT, dT = 0, 0, 0
all_proteins = ["CD28", "CD80", "CD86", "CTLA-4", "PD-1", "PD-L1", "PD-L2", "IL-2RA", "IL-12R",
                "IL-2", "IL-12"]

isValid = False
MCMC = 1
BD = 2


def show_proteins(window, check_buttons, boo_chosen_proteins):
    proteins_label = Label(window, text="Choose your proteins, and click 'Done':",
                           font=('Calibri bold', 18), background="lavender")
    proteins_label.place(x=300, y=280)
    for i in range(len(all_proteins)):
        c = Checkbutton(window, text=all_proteins[i], var=boo_chosen_proteins[i],
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
    iter_num.place(x=340, y=100)

    kT_label = Label(window, text="Enter kT value: ", font=('Calibri bold', 16),
                     background="lavender")
    kT_label.place(x=100, y=160)
    k = Entry(window, bd=2, font=('Calibri', 14), width=6)
    k.place(x=240, y=160)

    dT_label = Label(window, text="Enter Δt value: ", font=('Calibri bold', 16),
                     background="lavender")
    dT_label.place(x=100, y=220)
    d = Entry(window, bd=2, font=('Calibri', 14), width=6)
    d.place(x=240, y=220)
    return iter_num, k, d


def init_choices(choices_list):
    for i in range(len(all_proteins)):
        x = BooleanVar()
        x.set(False)
        choices_list.append(x)


def init_window(window, title, text, x, y):
    window.configure(background="lavender")
    window.geometry("1000x600")
    window.title(title)
    label = Label(window, text=text, font=("Calibri bold", 30),
                  background="lavender")
    label.place(x=x, y=y)


def BDgui(chosen_proteins):
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
            kT = int(k.get())
            dT = int(d.get())
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
    file_but.place(x=557, y=98)
    label = Label(window, text="Enter Δt value: ", font=('Calibri bold', 18),
                  background="lavender")
    d = Entry(window, bd=2, font=('Calibri', 14), width=6)
    label.place(x=300, y=200)
    d.place(x=460, y=205)
    choices_list = []
    init_choices(choices_list)
    check_buttons = []
    show_proteins(window, check_buttons, choices_list)

    def is_done():
        global dT, isValid
        try:
            dt_error, file_err = "Invalid Δt value", "Invalid file: must be csv file"
            dT = int(d.get())
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
    # sns.set_style('darkgrid')
    # sns.distplot([1,2,4,5,6])
    # plt.show()
    # Generate data on commute times.
    size, scale = 1000, 10
    commutes = pd.Series(np.random.gamma(scale, size=size) ** 1.5)
    print(commutes)
    commutes.plot.hist(grid=True, bins=20, rwidth=0.9,
                       color='#607c8e')
    # x = pd.Series(x)
    # x.plot.hist(grid=True, bins=20, rwidth=0.9, color='#607c8e')
    plt.title('Commute Times for 1,000 Commuters')
    plt.xlabel('Δt')
    plt.ylabel('Commute Time')
    plt.grid(axis='y', alpha=0.75)
    plt.show()


def main():
    global isValid
    algo = MAINgui()
    if algo == MCMC:
        run = MCMCgui
    else:
        run = BD
    while not isValid:
        chosen_proteins = []
        file_path = run(chosen_proteins)
    if algo == MCMC:
        return
    else:
        BD
    # histogram(0)


main()

from tkinter import *


def donothing():
    filewin = Toplevel(root)
    button = Button(filewin, text="Do nothing button")
    button.pack()


root = Tk()

root.title("SO - UFMG")
root.geometry(
    "500x250+200+200"
)

photo = PhotoImage(
    file="/home/rodrigo/Documentos/oitavoperiodo/optical_public/GUI/images/optical-fiber.png")
root.tk.call('wm', 'iconphoto',
                     root._w, photo)

menubar = Menu(root)
filemenu = Menu(menubar, tearoff=0)
filemenu.add_command(label="Novo", command=donothing)
filemenu.add_command(label="Abrir", command=donothing)
filemenu.add_command(label="Salvar", command=donothing)
filemenu.add_command(label="Salvar como...", command=donothing)
filemenu.add_command(label="Fechar", command=donothing)

filemenu.add_separator()

filemenu.add_command(label="Sair", command=root.quit)
menubar.add_cascade(label="Arquivo", menu=filemenu)
editmenu = Menu(menubar, tearoff=0)

editmenu.add_separator()

editmenu.add_command(label="Laser", command=donothing)
editmenu.add_command(label="Gerador de pulsos", command=donothing)
editmenu.add_command(label="Modulador", command=donothing)
editmenu.add_command(label="Fibra", command=donothing)
editmenu.add_command(label="Demodulador", command=donothing)

menubar.add_cascade(label="Componentes", menu=editmenu)
helpmenu = Menu(menubar, tearoff=0)
helpmenu.add_command(label="Help Index", command=donothing)
helpmenu.add_command(label="About...", command=donothing)
menubar.add_cascade(label="Help", menu=helpmenu)

root.config(menu=menubar)
root.mainloop()
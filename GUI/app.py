"""
Minha primeira GUI
"""
from re import M
from tkinter import *
from tkinter.ttk import *


menu_inicial = Tk()
menu_inicial.title("Primeiro app")
menu_inicial.geometry(
    "500x250+200+200"
)

menu_inicial.resizable(1,1)
# menu_inicial.minsize(
#     width=750,
#     height=750
# )

# menu_inicial.attributes('-zoomed', True)

import os 
print(os.getcwd())


photo = PhotoImage(
    file="/home/mfrodrigo/Documentos/ic/optical_public/GUI/images/fiber.png")
menu_inicial.tk.call('wm', 'iconphoto',
                     menu_inicial._w, photo)

menu_inicial['bg'] = 'White'

def cmd_click():
    print('Olá mundo !')
    
# botao 
cmd = Button(menu_inicial,
             text='Executar',
             command=cmd_click)
cmd.pack()

# resolução do nosso sistema 
largura_screen = menu_inicial.winfo_screenwidth()
altura_screen = menu_inicial.winfo_screenheight()

# posicao da janela
largura = 500
altura = 300
# posx = largura_screen/2 - largura/2
# posy = altura_screen/2 - altura/2
posx = 200
posy = 300

menu_inicial.geometry(f'{largura}x{altura}+{int(posx)}+{int(posy)}') 

label = Label(menu_inicial, text = "Label 1",
              font='Arial')
label.pack()
    
menu_inicial.mainloop()

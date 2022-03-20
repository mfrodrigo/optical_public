
# importing only those functions which
# are needed
from tkinter import Tk, mainloop, TOP
from tkinter.ttk import Button
  
# creating tkinter window
root = Tk()
  
# creating fixed geometry of the
# tkinter window with dimensions 150x200
root.geometry('200x150 + 400 + 300')
  
# Create Button and add some text
button = Button(root, text = 'Geeks')
button.pack(side = TOP, pady = 5)
  
# Execute Tkinter
root.mainloop()
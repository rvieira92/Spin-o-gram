#!/bin/env python
import tkinter as tk
from tkinter import ttk
from tkinter import colorchooser, filedialog
import argparse
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import glob

ROWS_DISP=10
FIG_SIZE=[5.4, 2.8]

# Auilixary functions (math+plot)
def rgba_to_hex(rgba):
    """Convert an RGBA tuple to a hexadecimal color string."""
    r, g, b, _ = rgba  # Ignore the alpha channel
    return "#{:02x}{:02x}{:02x}".format(int(r * 255),\
                                         int(g * 255),\
                                              int(b * 255))

def adjust_axis(ax):
    ax[0].set_xlabel(r'$\phi$')
    # Remove  labels from radial ticks
    ax[0].set_yticklabels([])
    # Sets the Zero on the top
    ax[0].set_theta_direction(-1)
    # Changes to clockwise
    ax[0].set_theta_offset(np.pi / 2.0)
    
    ax[0].set_thetalim(0,2*np.pi)
    ax[0].set_xticks(np.linspace(0,7*np.pi/4,8))
    ax[0].set_title("m$_{xy}$",loc='right')


    ax[1].set_xlabel(r'$\theta$')

    # Remove  labels from radial ticks
    ax[1].set_yticklabels([])
    
    ax[1].set_xticks(np.linspace(0,np.pi,5))
    # Sets the Zero on the top
    ax[1].set_theta_direction(-1)
    # Changes to clockwise
    ax[1].set_theta_offset(np.pi / 2.0)
    # Limits angles
    ax[1].set_thetalim(0,np.pi)
    ax[1].set_title("m$_z$",loc='right')

def get_momdist(restart_file,ntypes,NDIV=180) :

    raw_data = np.genfromtxt(restart_file) # cartesian
    iterens= np.unique(raw_data[:,0])[-1]
    iatom= np.unique(raw_data[:,1])[-1]
    idx = np.logical_and(raw_data[:,0]==iterens,raw_data[:,1]==iatom)
    #idx = np.logical_and(raw_data[:,0]==-1,raw_data[:,1]==1)
    raw_data = raw_data[idx][:,[3,4,5,6]] #

    # Convert to spherical coordinates according to the physical convention
    spherical = np.zeros((len(raw_data),3))    # Mag Mom magnitude
    spherical[:,0] = raw_data[:,0]    # phi, x-y plane
    spherical[:,1] = np.arctan2(raw_data[:,2], raw_data[:,1]) # theta elevation angle
    spherical[:,2] = np.arccos(raw_data[:,3]) 
    
    # Jacobian determinant to conserve density in non-linear transformation mz -> theta
    theta_weights =np.abs(1 / np.sqrt(1 - raw_data[:,3]**2)  )

    # Find different types of magnetic sites
    
    if ntypes == 0:
        mag_types = np.unique(spherical[:,0])
        ntypes = len(mag_types)
    
        # Get indexes of entries associated with each magnetic type
        idx_mag=[]
        for ii in range(ntypes):
            idx_mag.append(np.where(spherical[:,0]==mag_types[ii])[0])
    else :
        mag_types = np.arange(1,ntypes+1)

        idx_mag=[]
        for ii in range(ntypes):
            idx_mag.append(np.arange(ii,len(spherical),ntypes))


    cmap = plt.get_cmap('tab10') 
    colors = cmap(range(ntypes))
   
        
    phi_bins = np.linspace(-np.pi,np.pi,NDIV)
    theta_bins = np.linspace(0,np.pi,int(NDIV/2))
                                          
    # Initialize projected histogram arrays
    phi_hist_proj=np.zeros((ntypes,NDIV-1))
    theta_hist_proj=np.zeros((ntypes,int(NDIV/2)-1))

   # Count entries whithin specified bin intervals 
    for ii in range(ntypes):
        phi_hist_proj[ii]=np.histogram(spherical[idx_mag[ii],1],phi_bins)[0]
        phi_hist_proj[ii]/=phi_hist_proj[ii].max()
        theta_hist_proj[ii]=np.histogram(spherical[idx_mag[ii],2],\
                                         theta_bins,\
                                            weights=theta_weights[idx_mag[ii]])[0]
        theta_hist_proj[ii]/=theta_hist_proj[ii].max()


    # Get total histogram

    phi_hist_tot=np.histogram(spherical[:,1],phi_bins)[0]
    theta_hist_tot=np.histogram(spherical[:,2],theta_bins,weights=theta_weights)[0]

    phi_x =  np.concatenate(([phi_bins[0]],np.repeat(phi_bins[1:-1],2),[phi_bins[-1]]))
    theta_x = np.concatenate(([theta_bins[0]],np.repeat(theta_bins[1:-1],2),[theta_bins[-1]]))
 
    
    hist = [ {"phi": np.repeat(phi_hist_tot,2),"theta":np.repeat(theta_hist_tot,2),"color":colors[0]} ] 

    for ii in range(ntypes):
        hist.append({"label": mag_types[ii], "phi": np.repeat(phi_hist_proj[ii],2),\
                      "theta":np.repeat(theta_hist_proj[ii],2) ,"color": colors[ii]})
            
    bins={"phi":phi_x,"theta":theta_x}
  
    return bins,hist

# Auilixary functions (GUI)

def update_legend():
    """Update the legend to show only visible lines."""
    visible_lines = [line['theta'] for line in lines if line['theta'].get_visible()]
    # Update legend
    ax_proj[1].legend(handles=visible_lines,\
                      bbox_to_anchor=(1.04, 0.5),\
                          loc="center left",\
                              borderaxespad=0)
    canvas_proj.draw()

def toggle_visibility(line, var,canvas):
    """Toggle line visibility based on checkbox state."""
    value=var.get()
    line['phi'].set_visible(value)
    line['theta'].set_visible(value)
    canvas.draw()
    update_legend()

def pick_color(line, color_button,canvas):
    """Open a color picker and change line color."""
    color = colorchooser.askcolor(color=rgba_to_hex(line['phi'].get_facecolor()[0]))[1]
    
    if color:
        line['phi'].set_facecolor(color)
        line['phi'].set_color(color)
        line['theta'].set_facecolor(color)
        line['theta'].set_color(color)
        # Update the color of the associated button to match the new line color
        color_button.config(bg=color)
        canvas.draw()
        update_legend()

def save_plot(fig):
    """Save the plot to a file."""
    filepath = filedialog.asksaveasfilename(defaultextension=".png",\
                                            filetypes=[("PNG Files", "*.png"),\
                                                       ("PDF Files", "*.pdf"),\
                                                        ("All Files", "*.*")])
    if filepath:
        fig.savefig(filepath, bbox_inches="tight")
        print(f"Plot saved to {filepath}")

def rename_label(line, entry_widget):
    """Rename the line label and update the legend."""
    new_label = entry_widget.get()
    line['theta'].set_label(new_label)  # Update line label
    update_legend()  # Update the legend with the new label


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plots radial histogram of spins distribution')
    parser.add_argument('--file','-f', type=str,default='restart.*.out', help='UppASD output restart.*.out  file')
    parser.add_argument('--nsites','-s', type=int,default='0', help='Number of sites in posfile, changes proj mode')
    parser.add_argument('--NDIV', type=int, default='180', help='Number of bins to divide the radial angle in')
    args = parser.parse_args()
    
    restart_file = glob.glob(args.file)[0]
    nsites = args.nsites
    data_bins,data_hist= get_momdist(restart_file,nsites,args.NDIV)

    # Create the main window
    root = tk.Tk()
    root.title("Spin-o-gram")
    root.minsize(850, 450)
    #root.resizable(False, False)


    # Create notebook
    notebook = ttk.Notebook(root)
    notebook.pack(side=tk.TOP,expand=True, fill="both")
    # Create frames for notebook
    tab_tot = ttk.Frame(notebook)
    tab_proj = ttk.Frame(notebook)

    tab_proj.pack(fill="both",expand=True)
    tab_tot.pack(fill="both",expand=True)

    # Add frames to notebook
    notebook.add(tab_tot, text='Total')
    notebook.add(tab_proj, text='Projected')
  

    #############
    # Total Tab #
    #############

    # Create the Matplotlib figure and axis
    fig_tot, ax_tot = plt.subplots(nrows=1,ncols=2,subplot_kw={'polar':True},\
                                figsize=FIG_SIZE,dpi=125)
    
    line={'phi': ax_tot[0].fill_between(data_bins['phi'],\
                                        data_hist[0]['phi'],\
                                            color=data_hist[0]['color'],\
                                                alpha=0.65),\
          'theta': ax_tot[1].fill_between(data_bins['theta'],\
                                          data_hist[0]['theta'],\
                                            color=data_hist[0]['color'],\
                                                alpha=0.65)}
    
    adjust_axis(ax_tot)
    
    fig_tot.set_tight_layout(True)

    # Create a frame for the plot
    plot_frame = tk.Frame(tab_tot)
    plot_frame.pack(side=tk.RIGHT,expand=True, fill="both")

    canvas_tot = FigureCanvasTkAgg(fig_tot, master=plot_frame)
    canvas_tot.get_tk_widget().pack(expand=True, fill="both")


    control_panel_tot = tk.Frame(tab_tot)
    control_panel_tot.pack(side=tk.RIGHT,expand=False, fill="y",padx=10,pady=20)

    # Change color of plot
    panel_lbl = tk.Label(control_panel_tot, text='Change color:')
    panel_lbl.grid(row=0, column=0)

    color_button = tk.Button(control_panel_tot,\
                             bg=rgba_to_hex(line['phi'].get_edgecolor()[0]),\
                                command=lambda l=line,\
                                    cb=None: pick_color(l, cb,canvas_tot),\
                                          width=3, height=1)
    color_button.grid(row=1, column=0)

    # Update the lambda function to pass the button correctly
    color_button.config(command=lambda l=line, cb=color_button: pick_color(l, cb,canvas_tot))

    # Add a "Save Plot" button
    save_button_tot = tk.Button(control_panel_tot, text="Save Plot",command= lambda fig=fig_tot : save_plot(fig))
    save_button_tot.grid(row=2, column=0, pady=20)


    #################
    # Projected Tab #
    #################

    # Create a frame for the scrollbar(s).
    scroll_frame = tk.Frame(tab_proj)
    scroll_frame.pack(side=tk.LEFT,expand=False, fill="both")
    # Create labels for the columns
    panel_lbl = ttk.Label(scroll_frame,text='Show')
    panel_lbl.grid(row=0,column=0,ipady=10,sticky=tk.W)
    panel_lbl = ttk.Label(scroll_frame,text='Label')
    panel_lbl.grid(row=0,column=1,columnspan=2,ipady=10)
    panel_lbl = ttk.Label(scroll_frame,text='Color')
    panel_lbl.grid(row=0,column=3,ipady=10,sticky=tk.E)

    # Add a canvas in that frame.
    scroll_canvas = tk.Canvas(scroll_frame)
    scroll_canvas.grid(row=1, column=0,columnspan=4)

    # Create a vertical scrollbar linked to the canvas.
    vsbar = tk.Scrollbar(scroll_frame, orient=tk.VERTICAL, command=scroll_canvas.yview)
    vsbar.grid(row=1, column=4, sticky=tk.NS)
    scroll_canvas.configure(yscrollcommand=vsbar.set)

    # Create a frame on the canvas to contain the grid of buttons.
    projections_frame = tk.Frame(scroll_canvas)

    lines=[]
    # Create the Matplotlib figure and axis
    fig_proj, ax_proj = plt.subplots(nrows=1,ncols=2,subplot_kw={'polar':True},\
                                     figsize=FIG_SIZE,dpi=125)
    for ii in range(1,len(data_hist)):
        line = {'phi' : ax_proj[0].fill_between(data_bins['phi'],\
                                                data_hist[ii]['phi'],\
                                                    color=data_hist[ii]['color'],\
                                                        alpha=0.65),\
                'theta' : ax_proj[1].fill_between(data_bins['theta'],\
                                                  data_hist[ii]['theta'],\
                                                    color=data_hist[ii]['color'],\
                                                        label=data_hist[ii]['label'],\
                                                            alpha=0.65)}
        lines.append(line)
   

    ax_proj[1].legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
    adjust_axis(ax_proj)
    fig_proj.set_tight_layout(True)
    for i, line in enumerate(lines):
        
        var = tk.BooleanVar(value=True)
        checkbox = tk.Checkbutton(projections_frame, variable=var, command=lambda l=line, v=var: toggle_visibility(l, v,canvas_proj))
        checkbox.grid(row=i+1, column=0) 
            
        # Entry to rename the label
        label_entry = tk.Entry(projections_frame)
        label_entry.insert(0, line['theta'].get_label())
        label_entry.grid(row=i+1, column=1)
        
        # Button to rename label
        rename_button = tk.Button(projections_frame, text="Rename", command=lambda l=line, entry=label_entry: rename_label(l, entry))
        rename_button.grid(row=i+1, column=2)

        color_button = tk.Button(projections_frame,bg=rgba_to_hex(line['phi'].get_edgecolor()[0]), command=lambda l=line, cb=None: pick_color(l, cb,canvas_proj), width=2, height=1)
        color_button.grid(row=i+1, column=3,padx=10)
        color_button.config(command=lambda l=line, cb=color_button: pick_color(l, cb,canvas_proj))

    # Create canvas window to hold the buttons_frame.
    scroll_canvas.create_window((0,0), window=projections_frame, anchor=tk.NW)

    projections_frame.update_idletasks()  # Needed to make bbox info available.
    bbox = scroll_canvas.bbox(tk.ALL)  # Get bounding box of canvas with Buttons.

    # Define the scrollable region as entire canvas with only the desired
    # number of rows and columns displayed.
    w, h = bbox[2]-bbox[1], bbox[3]-bbox[1]
    ROWS=len(lines)
    dh = int((h/ROWS) * ROWS_DISP)
    scroll_canvas.configure(scrollregion=bbox, height=dh,width=w)
    scroll_canvas.configure(scrollregion=scroll_canvas.bbox("all"))

    # Add a "Save Plot" button
    save_button_proj = tk.Button(scroll_frame, text="Save Plot",command= lambda fig=fig_proj : save_plot(fig))
    save_button_proj.grid(row=2, column=0,columnspan=4, pady=10)

    # Create a frame for the plot
    plot_frame = tk.Frame(tab_proj)
    plot_frame.pack(side=tk.RIGHT,expand=True, fill="both")


    canvas_proj = FigureCanvasTkAgg(fig_proj, master=plot_frame)
    canvas_proj.get_tk_widget().pack(expand=True, fill="both")

    #################
    #################

    #Add footer
    footer = ttk.Frame(root)
    footer.pack(side=tk.TOP,expand=False, fill="both")
    footer_txt = ttk.Label(footer,anchor="e",justify=tk.RIGHT,\
                        text='by Rafael Vieira\nUppsala University')
    footer_txt.pack(side='right',pady=10,padx=10)
    footer_txt = ttk.Label(footer,anchor="e",justify=tk.RIGHT,\
                        text='[Spin-o-gram] : Spin distribution plotter for UppASD')
    footer_txt.pack(side='left',pady=10,padx=10)

    # Start the Tkinter main loop
    root.mainloop()



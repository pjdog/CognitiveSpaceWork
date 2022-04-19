# -*- coding: utf-8 -*-
"""
Orbital Dynamics Homework Problem
Created on Sat Apr 16 13:19:40 2022

@author: Paul Hughes
"""

"""
GUI Template for 
Keplerian User Inputs:
"""
from tkinter import * 
from tkinter.ttk import Progressbar
from ttkthemes import  ThemedTk
import pandas as pd
import numpy as np
import os
import sys
import psutil
import logging
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import PIL
import load_constants as CON

#from pylab import *

"""
Button Behavior
"""
def submit():
    a_st = a_in.get()
    e_st = e_in.get()
    i_st = i_in.get()
    w_st = w_in.get()
    O_st = O_in.get()
    v_st = v_in.get()
    t_st = t_in.get()
    inputChecks(a_st,e_st,i_st,w_st,O_st,v_st,t_st)

def delete():
    a_in.delete(0,END)
    e_in.delete(0,END)
    i_in.delete(0,END)
    w_in.delete(0,END)
    O_in.delete(0,END)
    v_in.delete(0,END)
    t_in.delete(0,END)
    focus_out_a("")
    focus_out_i("")
    focus_out_e("")
    focus_out_O("")
    focus_out_w("")
    focus_out_v("") 
    focus_out_t("")
    
def degI_toggle():
    if degInc_btn.config('relief')[-1] == 'sunken':
        degInc_btn.config(relief="raised")
        degInc_btn.config(text = "Degrees")
    else:
        degInc_btn.config(relief="sunken")
        degInc_btn.config(text = "Radians")

    
def degw_toggle():
    if degw_btn.config('relief')[-1] == 'sunken':
        degw_btn.config(relief="raised")
        degw_btn.config(text = "Degrees")
    else:
        degw_btn.config(relief="sunken")
        degw_btn.config(text = "Radians")
        
            
def degO_toggle():
    if degO_btn.config('relief')[-1] == 'sunken':
        degO_btn.config(relief="raised")
        degO_btn.config(text = "Degrees")
    else:
        degO_btn.config(relief="sunken")
        degO_btn.config(text = "Radians")
        
            
def degv_toggle():
    if degv_btn.config('relief')[-1] == 'sunken':
        degv_btn.config(relief="raised")
        degv_btn.config(text = "Degrees")
    else:
        degv_btn.config(relief="sunken")
        degv_btn.config(text = "Radians")

def anomaly_type():
    if anomaly_btn.config('text')[-1]=='True Anomaly':
        anomaly_btn.config(text = 'Mean Anomaly')
    elif anomaly_btn.config('text')[-1]=='Mean Anomaly':
        anomaly_btn.config(text = 'Eccentric Anomaly')
    else:
        anomaly_btn.config(text = 'True Anomaly')

def time_span():
    if time_btn.config('text')[-1]=='Days':
        time_btn.config(text = 'Hours')
    elif time_btn.config('text')[-1]=='Hours':
        time_btn.config(text = 'Minutes')
    elif time_btn.config('text')[-1]=='Minutes':
        time_btn.config(text = 'Seconds')
    else:
        time_btn.config(text = 'Days')
'''
Handling Default Entries
'''
def focus_in(btn,string):
    if btn.get() == string:
        btn.delete(0,END)
        btn.config(fg='black')
def focus_in_a(_):
    if a_in.get() == "6730.773":
        a_in.delete(0,END)
        a_in.config(fg='black')
def focus_out_a(_):
    if a_in.get() == "":
        a_in.delete(0,END)
        a_in.config(fg = 'grey')
        a_in.insert(0,"6730.773")

def focus_in_e(_):
    if e_in.get() == "0.0008685":
        e_in.delete(0,END)
        e_in.config(fg='black')
def focus_out_e(_):
    if e_in.get() == "":
        e_in.delete(0,END)
        e_in.config(fg = 'grey')
        e_in.insert(0,"0.0008685")

def focus_in_i(_):
    if i_in.get() == "51.6408":
        i_in.delete(0,END)
        i_in.config(fg='black')
def focus_out_i(_):
    if i_in.get() == "":
        i_in.delete(0,END)
        i_in.config(fg = 'grey')
        i_in.insert(0, "51.6408")


def focus_in_w(_):
    if w_in.get() == "246.9642":
        w_in.delete(0,END)
        w_in.config(fg='black')
def focus_out_w(_):
    if w_in.get() == "":
        w_in.delete(0,END)
        w_in.config(fg = 'grey')
        w_in.insert(0,"246.9642")



def focus_in_O(_):
    if O_in.get() =="261.4942":
        O_in.delete(0,END)
        O_in.config(fg='black')
def focus_out_O(_):
    if O_in.get() =="":
        O_in.delete(0,END)
        O_in.config(fg = 'grey')
        O_in.insert(0, "261.4942")


def focus_in_v(_):
    if v_in.get() == "257.7333":
        v_in.delete(0,END)
        v_in.config(fg='black')
def focus_out_v(_):
    if v_in.get() == "":
          v_in.delete(0,END)
          v_in.config(fg = 'grey')
          v_in.insert(0,"257.7333")
    
def focus_in_t(_):
    if t_in.get() == "1":
        t_in.delete(0,END)
        t_in.config(fg='black')
def focus_out_t(_):
    if t_in.get() == "":
          t_in.delete(0,END)
          t_in.config(fg = 'grey')
          t_in.insert(0,"1")
    

def restart_sim():
    """Restarts the sim in case of invalid input, with file objects and descriptors
       cleanup
    """

    try:
        p = psutil.Process(os.getpid())
        for handler in p.get_open_files() + p.connections():
            os.close(handler.fd)
    except Exception as e:
        logging.error(e)

    python = sys.executable
    os.execl(python, python, *sys.argv)

'''
Math Scripts
'''
def inputChecks(a_st,e_st,i_st,w_st,O_st,v_st, t_st):
    ''' 
    This function takes the inputs and makes sure they're in the right format and then changes them to 
    radians, and true anomaly in order to then change to ECEF Coordinates
    '''
    try:
        a = float(a_st)
        

    except:
        popup = Tk()
        popup.wm_title("!")
        warn_label =Label(popup, text = "Error: Unacceptable Character in Semimajor Axis")
        warn_label.pack(side = 'top')
        restart_btn = Button(popup, text = "Restart Simulation",command = restart_sim)
        restart_btn.pack(side = 'bottom')
    
    try:
        e = float(e_st)
        if 1<=e:
            popup = Tk()
            popup.wm_title("!")
            warn_label =Label(popup, text = "Unfortunately parabolic and hyperbolic orbits have not been considered")
            warn_label.pack(side = 'top')
            restart_btn = Button(popup, text = "Restart Simulation",command = restart_sim)
            restart_btn.pack(side = 'bottom')

    except:
        popup = Tk()
        popup.wm_title("!")
        warn_label =Label(popup, text = "Error: Unacceptable Character in Eccentricity")
        warn_label.pack(side = 'top')
        restart_btn = Button(popup, text = "Restart Simulation",command = restart_sim)
        restart_btn.pack(side = 'bottom')


    try:
        i = float(i_st)
        if degInc_btn.config('relief')[-1] == 'raised':
            i = m_deg2rad*i
        if i<-np.pi/2.0 or i>np.pi/2.0:
            popup = Tk()
            popup.wm_title("!")
            warn_label =Label(popup, text = "Error: Inclination must be between +/- 90 degrees")
            warn_label.pack(side = 'top')
            restart_btn = Button(popup, text = "Restart Simulation",command = restart_sim)
            restart_btn.pack(side = 'bottom')
    except:
        popup = Tk()
        popup.wm_title("!")
        warn_label =Label(popup, text = "Error: Unacceptable Character in inclination")
        warn_label.pack(side = 'top')
        restart_btn = Button(popup, text = "Restart Simulation",command = restart_sim)
        restart_btn.pack(side = 'bottom')
        

    try:
        w = float(w_st)
        if degw_btn.config('relief')[-1] == 'raised':
             w= m_deg2rad*w
    except:
        popup = Tk()
        popup.wm_title("!")
        warn_label =Label(popup, text = "Error: Unacceptable Character in Argument of Perigee")
        warn_label.pack(side = 'top')
        restart_btn = Button(popup, text = "Restart Simulation",command = restart_sim)
        restart_btn.pack(side = 'bottom')
    
    try:
        O = float(O_st)
        if degO_btn.config('relief')[-1] == 'raised':
             O= m_deg2rad*O

    except:
        popup = Tk()
        popup.wm_title("!")
        warn_label =Label(popup, text = "Error: Unacceptable Character in Right Ascension of Ascending Node")
        warn_label.pack(side = 'top')
        restart_btn = Button(popup, text = "Restart Simulation",command = restart_sim)
        restart_btn.pack(side = 'bottom')
        
    
    try:
        v = float(v_st)
        if degv_btn.config('relief')[-1] == 'raised':
             v= m_deg2rad*v
    except:
        popup = Tk()
        popup.wm_title("!")
        warn_label =Label(popup, text = "Error: Unacceptable Character in Anomaly Term")
        warn_label.pack(side = 'top')
        restart_btn = Button(popup, text = "Restart Simulation",command = restart_sim)
        restart_btn.pack(side = 'bottom')
    if anomaly_btn.config('text')[-1]=='Eccentric Anomaly':
        v = eccentric2true(v,e)
    elif   anomaly_btn.config('text')[-1]=='Mean Anomaly':
        v = mean2true(v,e)
    v = abs(v)
    try:
        tf = float(t_st)
        if tf<0:
            popup = Tk()
            popup.wm_title("!")
            warn_label =Label(popup, text = "Error: Negative Time")
            warn_label.pack(side = 'top')
            restart_btn = Button(popup, text = "Restart Simulation",command = restart_sim)
            restart_btn.pack(side = 'bottom')

        if time_btn.config('text')[-1]=='Days':
            tf=3600*tf*24
        elif time_btn.config('text')[-1]=='Hours':
            tf=tf*3600*tf
        elif time_btn.config('text')[-1]=='Minutes':
            tf =tf*60
    except:
        popup = Tk()
        popup.wm_title("!")
        warn_label =Label(popup, text = "Error: Unacceptable Character in Anomaly Term")
        warn_label.pack(side = 'top')
        restart_btn = Button(popup, text = "Restart Simulation",command = restart_sim)
        restart_btn.pack(side = 'bottom')
    '''
    Initialize Simulation:
        Now that the input data has been verified the simulation will be initialized
    '''
    cart0 = np.zeros(6)
    cart0 = kep2cart(a,e,i,O,w,v);
    
    window.destroy()
    loading =Tk()
    loading.wm_title("Simulating!")
    progress_var = DoubleVar()
    progress= Progressbar(loading,variable = progress_var,maximum = tf, orient = HORIZONTAL, length = 300, mode = 'determinate')

    top_label =Label(loading, text = "Simulation in progress")
    top_label.pack(side = 'top')
    # initialize data frame for storing info for quicker writing to file or manipulation later
    progress.pack(pady=10)
    cart1=cart0
    rn =np.linalg.norm(np.array([cart0[0],cart0[1],cart0[2]]))
    vn = np.linalg.norm(np.array([cart0[3],cart0[4],cart0[5]]))
    sim_out = pd.DataFrame(
        data = {'Time':0.0,
                'x': cart0[0],
                'y': cart0[1],
                'z': cart0[2],
                'vx': cart0[3],
                'vy': cart0[4],
                'vz': cart0[5],
                'rn': rn,
                'vn': vn,
                'a': a,
                'e': e,
                'i': i,
                'w': w,
                'O': O,
                'v': v
                },
    index = [0],
    columns  =['Time','x','y','z','vx','vy','vz','rn','vn','a','e','i','w','O','v'])
    for t_k in np.arange(0.0,tf,CON.delta_t)[1:]:
        cart1 = rk45(t_k,cart0)
        #get our keplerian elementrs back if we want to print
        if t_k% CON.save_data <CON.delta_t:
            #print('got to')
            #print(t_k)
            kep = cart2kep(cart1[:3].flatten(),cart1[3:7].flatten())
            point=  {'Time':t_k,
                'x': cart1 [0,0],
                'y': cart1[1,0],
                'z': cart1[2,0],
                'vx': cart1[3,0],
                'vy': cart1[4,0],
                'vz': cart1[5,0],
                'rn': np.linalg.norm(np.array([cart1[0,0],cart1[1,0],cart1[2,0]])),
                'vn': np.linalg.norm(np.array([cart1[3,0],cart1[4,0],cart1[5,0]])),
                'a': kep[0],
                'e': kep[1],
                'i': kep[2],
                'w': kep[4],
                'O': kep[3],
                'v': kep[5]
                }
            sim_out.loc[int(t_k)] = point
            progress_var.set(t_k)
            loading.update()
        cart0 = cart1
        
    
    #progress.stop
    loading.destroy()
    printer_dialogue(sim_out)

    
#Transforming between anomalies
def eccentric2true(E,e):
    return 2.0*np.arctan(np.sqrt((1.0+e)/(1.0-e))*np.tan(E/2.0))
def mean2true(M,e):
    error_tolerance = 1e-8;
    if M<np.pi:
        E=M + e/2;
    else:
        E = M-e/2;
    ratio = 1.0
    n = 0;
    #Using newtons method w/ max of 100 poss steps:
    while abs(ratio)>error_tolerance and n <100:
        ratio = (E-e*np.sin(E)-M)/(1-e*np.cos(E));
        E = E - ratio
        n = n+1
    return eccentric2true(E,e)
            
'''
State Vector Transformation
'''
def kep2cart(a,e,i,O,w,v):
    p = a*( 1.0-np.power(e,2))
    #Dealing with edge cases
    if abs(i)<0.01 and abs(e)<0.01:
        v = v + w+O
        w = 0.0
        O = 0.0
    if abs(e)<0.01 and abs(i)>0.01:
        v = v+w+O
        w =0.0
    rpqw = np.array([[p*np.cos(v)/(1+e*np.cos(v))],[p*np.sin(v)/(1+e*np.cos(v))],[0] ]      )
    vpqw = np.array([[-np.sqrt(CON.mu_E/p)*np.sin(v)],[np.sqrt(CON.mu_E/p)*(e+np.cos(v))],[0]])
    rot = np.matmul(rot3(-O),np.matmul(rot1(-i),rot3(-w)))
    rIJK = rot.dot(rpqw)
    vIJK = rot.dot(vpqw)
    return np.concatenate((rIJK,vIJK),axis = 0)
def rot3(x):
    return np.array([[np.cos(x), np.sin(x),0],
                     [-1.0*np.sin(x), np.cos(x),0],
                     [0,0,1.0]])
def rot1(x):
    return np.array([[1.0,0,0],
                     [0,np.cos(x),np.sin(x)],
                     [0,-1.0*np.sin(x),np.cos(x)]])     
def cart2kep(r,v):
    h = np.cross(r,v)
    hn = np.linalg.norm(h)
    n = np.cross(np.array([0,0,1.0]),h)
    
    nn = np.linalg.norm(n)
    n = n*(1/nn)
    rn = np.linalg.norm(r)
    vn = np.linalg.norm(v)
    e = ((np.power(vn,2)-CON.mu_E/rn)*r-(r.dot(v))*(v))/CON.mu_E
    en = np.linalg.norm(e)
    xsi = np.power(vn,2)/2.0-CON.mu_E/rn
    if en != 1:
        a = -CON.mu_E/(2.0*xsi)
        p = a*(1-np.power(e,2))
    else:
        p = np.pow(h,2)/CON.mu_E
    i = np.arccos(h[2]/hn)
    O = np.arccos(n[0]/nn)
    if n[1]<0:
        O = np.pi*2.0-O
    w = np.arccos(n.dot(e)/(en*nn))  
    if e[2]<0:
        w = np.pi*2.0-w
    v = np.arccos(e.dot(r)/(en*rn))
    return np.array([a,en,i,O,w,v])

    
def dataPlot(data, x_variable, y_variable, save):
    x_var= x_variable.get()
    y_var = y_variable.get()
    graphTitle = y_var+" as a function of "+ x_var
    plt_2D= data.plot(x= x_var,y= y_var, title = graphTitle)
    plt_2D.set_xlabel(x_var)
    plt_2D.set_ylabel(y_var)
    if save:
        fig = plt_2D.get_figure()
        fig.savefig(graphTitle+'.png')

def plot_earth(data, showEarth, showSpeedMap, az, el): 
        showEarth= showEarth.get()
        showSpeedMap = showSpeedMap.get()
        az = az.get()
        el = el.get()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        if showEarth:
            # load bluemarble with PIL
            bm = PIL.Image.open('earth_map.jpg')
            # make the image smaller
            bm = np.array(bm.resize([d//5 for d in bm.size]))/256
            lons = np.linspace(-180, 180, bm.shape[1]) * np.pi/180 
            lats = np.linspace(-90, 90, bm.shape[0])[::-1] * np.pi/180 
            x = np.outer(np.cos(lons), np.cos(lats)).T
            y = np.outer(np.sin(lons), np.cos(lats)).T
            z = np.outer(np.ones(np.size(lons)), np.sin(lats)).T
            #This isn't accurate to the obliquity of the earth but it's fun
            r_E = CON.r_E
            ax.plot_surface(x*r_E, y*r_E, z*r_E, rstride=4, cstride=4, facecolors = bm)

        x_data = data.loc[:,"x"].tolist()
        y_data = data.loc[:,"y"].tolist()
        z_data = data.loc[:,"z"].tolist()
        v_data = data.loc[:,'vn'].tolist()
        img = ax.scatter(x_data,y_data,z_data,c = v_data, cmap = 'hsv')
        if showSpeedMap:
            fig.colorbar(img,ax=ax)
        ax.view_init(az, el)
        fig.savefig("plottedSatellite.png")

def printer_dialogue(data):
    plot_dialogue = Tk()
    OPTIONS = data.columns.tolist()
    x_variable = StringVar()
    y_variable = StringVar()
    showEarth = BooleanVar()
    showSpeedMap = BooleanVar()
    az= DoubleVar()
    el = DoubleVar()
    showEarth.set(False)
    showSpeedMap.set(True)
    az.set(-60.0)
    el.set = (30.0)
    w_az = Scale(plot_dialogue, from_=-180, to=180, orient=HORIZONTAL,variable =az)
    w_el = Scale(plot_dialogue, from_=-90, to=90, orient=HORIZONTAL,variable = el)

    x_variable.set(OPTIONS[0])
    
    y_variable.set(OPTIONS[1])
    w_x = OptionMenu(plot_dialogue, x_variable,*OPTIONS)
    w_x.config(text=x_variable)
    w_y = OptionMenu(plot_dialogue, y_variable,*OPTIONS)
    w_y.config(text = y_variable)
    Label(plot_dialogue, text = "Graph 2D Variables").grid(row =1)
    w_x.grid(row = 1, column =1)
    w_y.grid(row= 1, column = 2)
    graph_btn = Button(plot_dialogue,text = "plot 2D Data", command = lambda :dataPlot(data, x_variable, y_variable,False))
    graph_btn.grid(row = 1, column = 3)
    save_graph_btn= Button(plot_dialogue,text = "save plot", command = lambda :dataPlot(data, x_variable, y_variable,True))
    save_graph_btn.grid(row = 1, column = 4)
    save_csv_btn = Button(plot_dialogue, text = "save csv of sim data", command= data.to_csv('out.csv', index = False))
    save_csv_btn.grid(row = 6, column = 2)
    Label(plot_dialogue, text = "Configure 3D Plot Output:").grid(row=2, columnspan = 4)
    earth_graph_btn = Button(plot_dialogue,text = "Export 3D Plot", command = lambda :plot_earth(data, showEarth, showSpeedMap, az, el))
    earth_graph_btn.grid(row = 6, column = 1)
    earth_btn= Checkbutton(text="Show Earth",variable = showEarth, onvalue = True, offvalue = False)
    showSpeedMap_btn= Checkbutton(text="Show Speed Map",variable = showSpeedMap, onvalue = True, offvalue = False)
    earth_btn.grid(row = 3, column =4)
    showSpeedMap_btn.grid(row = 4, column = 4)
    Label(plot_dialogue, text = "Azimuthal Viewing Angle").grid(row = 3)
    w_az.grid(row=3, column=1,columnspan = 2)
    Label(plot_dialogue, text = "Elevation Viewing Angle").grid(row = 4)
    w_el.grid(row=4, column=1,columnspan = 2)
'''
Acceleration vector
'''
def accel(t,cart):
    out = np.zeros((6,1))
    out[0,0] = cart[3,0]
    out[1,0] = cart[4,0]
    out[2,0] = cart[5,0]
    rn = np.linalg.norm(np.array([cart[0,0],cart[1,0],cart[2,0]]))
    #two body acceleration 
    a2x = -CON.mu_E/np.power(rn,3)*cart[0,0]
    a2y = -CON.mu_E/np.power(rn,3)*cart[1,0]
    a2z = -CON.mu_E/np.power(rn,3)*cart[2,0]
    if J2_switch== 1:
        a2x=a2x+-3*J2*np.power(r_E,2)*CON.mu_E*cart[0,0]/(2*rn^5)*(1-5*np.power(cart[2]))/np.power(rn,2);
        a2y=a2y+-3*J2*np.power(r_E,2)*CON.mu_E*cart[1,0]/(2*rn^5)*(1-5*np.power(cart[2]))/np.power(rn,2);
        a2z=a2z+-3*J2*np.power(r_E,2)*CON.mu_E*cart[2,0]/(2*rn^5)*(3-5*np.power(cart[2]))/np.power(rn,2);
    if J345_switch== 1:
    #J2 acceleration 
        a2x=a2x+0;
        a2y=a2y+0;
        a2z=a2z+0; 
    out[3,0] =a2x
    out[4,0] =a2y
    out[5,0] =a2z
    return out
    #return(np.concatenate(cart[35])
'''
Integration:
    A simple runge-kutta 4th orcdf
'''
def rk45(t,cart0):
    #f' 
    k1 = accel(t,cart0)
    k2 = accel(t+0.5*CON.delta_t,cart0+0.5*k1*CON.delta_t)
    k3 = accel(t+0.5*CON.delta_t,cart0+0.5*k2*CON.delta_t)
    k4 = accel(t+ CON.delta_t, cart0+k3*CON.delta_t)
    return cart0+ CON.delta_t/6.0*(k1+2*k2+2*k3+k4)
    

'''
Program Start:
'''  
#Constants
m_deg2rad = np.pi/180.0
#GUI STARTUP  
window = ThemedTk(theme='arc')
window.title('Run Config:')
'''
Button Setup:
'''

submit_btn = Button(window,text = "submit simulation", command = submit)
delete_btn = Button(window,text = "clear inputs",command = delete)
degInc_btn = Button(text="Degrees", width=12, relief="raised", command=degI_toggle)
degw_btn = Button(text="Degrees", width=12, relief="raised", command=degw_toggle)
degO_btn = Button(text="Degrees", width=12, relief="raised", command=degO_toggle)
degv_btn = Button(text="Degrees", width=12, relief="raised", command=degv_toggle)
time_btn = Button(text="Days",width = 12,command = time_span)
anomaly_btn = Button(text="Mean Anomaly",width = 20, command = anomaly_type)
J2_switch=IntVar()
J345_switch =IntVar()
J2_btn= Checkbutton(text="J2 Perturbations",variable = J2_switch, onvalue = 1, offvalue = 0)
J345_btn = Checkbutton(text ="J3, J4, and J5 Perturbations",variable = J345_switch, onvalue = 1, offvalue = 0)

'''
Style Setup
'''
#delete.pack(side = BOTTOM)
#entry = Entry()
#entry.config(font =( 'Arial', 25))
#entry.config(bg='black')
#entry.config(fg = 'white')
#entry.insert(0,'Please input your initial keplerian orbital elements:')
'''
Input Setup:
'''
Label(window, text="Please input initial Keplerian Elements, or press submit to simulate the ISS:").grid(row =0,columnspan = 3,sticky='EW')
Label(window, text="SemiMajor Axis:").grid(row=1)
Label(window, text="km").grid(row=1,column =2)

Label(window, text="Eccentricity:").grid(row=2)
Label(window, text="Inclination:").grid(row=3)
Label(window, text="Right Ascencion of the Node:").grid(row=4)
Label(window, text="Argument of Perigee:").grid(row=5)
Label(window, text="Time of propagation: ").grid(row =7)

'''
Organizational Setup:
'''
a_in = Entry(window, bg='white', width=30, fg='grey')
e_in = Entry(window, bg='white', width=30, fg='grey')
i_in = Entry(window, bg='white', width=30, fg='grey')
O_in = Entry(window, bg='white', width=30, fg='grey')
w_in = Entry(window, bg='white', width=30, fg='grey')
v_in = Entry(window, bg='white', width=30, fg='grey')
t_in = Entry(window, bg='white', width=30,fg='grey')
a_in.grid(row=1, column=1)
a_in.insert(0,"6730.773")
a_in.bind("<FocusIn>",focus_in_a)
a_in.bind("<FocusOut>",focus_out_a)
e_in.grid(row=2, column=1)
e_in.insert(0,"0.0008685")
e_in.bind("<FocusIn>",focus_in_e)
e_in.bind("<FocusOut>",focus_out_e)
i_in.grid(row=3, column =1)
i_in.insert(0, "51.6408")
i_in.bind("<FocusIn>",focus_in_i)
i_in.bind("<FocusOut>",focus_out_i)
O_in.grid(row=4, column =1)
O_in.insert(0, "261.4942")
O_in.bind("<FocusIn>",focus_in_O)
O_in.bind("<FocusOut>",focus_out_O)
w_in.grid(row=5, column =1)
w_in.insert(0,"246.9642")
w_in.bind("<FocusIn>",focus_in_w)
w_in.bind("<FocusOut>",focus_out_w)
v_in.grid(row=6, column =1)
v_in.insert(0,"257.7333")
v_in.bind("<FocusIn>",focus_in_v)
v_in.bind("<FocusOut>",focus_out_v)
t_in.grid(row=7,column = 1)
t_in.insert(0,'1')
t_in.bind("<FocusIn>",focus_in_t)
t_in.bind("<FocusOut>",focus_out_t)

degInc_btn.grid(row=3,column =2)
degO_btn.grid(row =4,column =2) 
degw_btn.grid(row =5,column =2) 
degv_btn.grid(row =6,column =2) 
anomaly_btn.grid(row =6, column = 0)
time_btn.grid(row =7, column =2)
J2_btn.grid(row=8, column = 1)
J345_btn.grid(row=9, column= 1)
delete_btn.grid(row=10, column = 0)
submit_btn.grid(row =10, column = 2)

window.mainloop()


#----------------------------------------------------------------------------------------------------
#----------------------------------------LIBRARIES FOR PROJECT---------------------------------------
#----------------------------------------------------------------------------------------------------
import pyaudio
import struct
import math
import wave
import random
from scipy import signal
from matplotlib import pyplot
from myfunctions import clip16
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import tkinter as Tk

# Define Tkinter root
root = Tk.Tk()
root.geometry("1200x500")
root.title('PHASE VOCODER BY MUNIBA & JEEVTESH')

#----------------------------------------------------------------------------------------------------
#----------------------------------------FUNCTIONS FOR PROJECT---------------------------------------
#----------------------------------------------------------------------------------------------------
def CmplxAng(ang):
    pphi = np.ndarray(shape=(1,1), dtype=complex)
    phi = math.cos(ang)+1j*math.sin(ang)
    return(phi)

def numGen(lenx, leny):
    phi = np.ndarray(shape=(lenx,leny), dtype=complex)
    G = 1
    for i in range(leny):
                     for j in range(lenx):
                         endR = int(math.radians(360))
                         rand_ang = random.randrange(0,endR)
                         phi[j,i] = G*(math.cos(rand_ang)+1j*math.sin(rand_ang))
    return(phi)

def Normal():
  global a, b, c, d, e
  c = True
  b = False
  a = False
  d = False
  e = False

def Robot():
  global a, b, c, d, e
  c = False
  b = True
  a = False
  d = False
  e = False

def Whisper():
    global a, b, c, d, e
    c = False
    b = False
    a = True
    d = False
    e = False

def TimeSt():
    global a, b, c, d, e
    c = False
    b = False
    a = False
    d = True
    e = False

# Update plot when slider is moved	
def hopUpdate(event):
    global ALen
    global SLen
    global hPr
    hPr = HOP.get()
    SLen = hPr * ALen
    ALen = SLen / hPr

def DeNoise():
    global a, b, c, d, e
    c = False
    b = False
    a = False
    d = False
    e = True

def VolUpdate(event):
    global gain
    gain = VOL.get()

def Quit():
  global QUIT
  print('Good Bye')
  QUIT = True
  root.quit()

a = False
b = False
c = True    #NORMAL
d = False
e = False
QUIT = False

#----------------------------------------------------------------------------------------------------
#-------------------------------VARIABLES FOR BLOCK LENGTH AND STFT----------------------------------
#----------------------------------------------------------------------------------------------------
#for real time on computer speakers BLOCKLEN = 4096, RATE = 44100
#for real time on headphones BLOCKLEN = 1024, RATE = 32000

BLOCKLEN = 1024      # Number of frames per block
WIDTH = 2           # Number of bytes per signal value
CHANNELS = 1        # mono
RATE = 16000        # Frame rate (frames/second)
RECORD_SECONDS = 5000

stft_N = BLOCKLEN
R_lap = 0
MAXVALUE = 2**15-1
gain = 1
zxx_len = int((stft_N/2)+1)

#Pitch Shifting
ALen = 80
SLen = 80

hPr = SLen/ALen

flag = 0
loop1 = True
#----------------------------------------------------------------------------------------------------
#----------------------------------------PYAUDIO-----------------------------------------------------
#----------------------------------------------------------------------------------------------------
# Open the audio output stream
p = pyaudio.PyAudio()

stream = p.open(
    format      = p.get_format_from_width(WIDTH),
    channels    = CHANNELS,
    rate        = RATE,
    input       = True,
    output      = True)
num_blocks = int(RATE / BLOCKLEN * RECORD_SECONDS)
#----------------------------------------------------------------------------------------------------
#----------------------------------------PYPLOT------------------------------------------------------
#----------------------------------------------------------------------------------------------------
pyplot.ion()           # Turn on interactive mode so plot gets updated

output_block = [0] * BLOCKLEN
t = [n*1000/float(RATE) for n in range(BLOCKLEN)]


fig, ax = pyplot.subplots()
im = pyplot.imread("background2.png")
my_plot = fig.add_subplot(1, 1, 1)
ax.imshow(im)
ax.axis('off')
[g2] = my_plot.plot(t, output_block, color = '#211A81', linewidth = 0.6)
#blue color = #211A81, gray color = 505051



my_plot.set_xlim(0, 1000.0 * BLOCKLEN/RATE)         # set x-axis limits
my_plot.set_ylim(-3500, 3500)        # set y-axis limits
my_plot.set_xlabel('Time (milliseconds)')
my_plot.axis('off')


#----------------------------------------------------------------------------------------------------
#--------------------------------------------GUI INFO------------------------------------------------
#----------------------------------------------------------------------------------------------------

# Define Tk variables
HOP = Tk.DoubleVar()
VOL = Tk.DoubleVar()

# Define widgets
B_n = Tk.Button(root, text = 'Normal', command = Normal, font=("Copperplate"), fg = 'white',highlightbackground='#0059b3')
B_r = Tk.Button(root, text = 'Robotization', command = Robot, font=("Copperplate"), fg = 'white',highlightbackground='#0059b3')
B_w = Tk.Button(root, text = 'Whisperization', command = Whisper, font=("Copperplate"), fg = 'white',highlightbackground='#0059b3')
B_t = Tk.Button(root, text = 'Pitch Shifting', command = TimeSt, font=("Copperplate"), fg = 'white',highlightbackground='#0059b3')
B_quit = Tk.Button(root, text = 'Quit', command = Quit, font=("Copperplate"), fg = 'white',highlightbackground='#0059b3')
L1 = Tk.Label(root, text = 'PHASE VOCODER', font=("Lucida Grande", 38), fg = '#0059b3')
B_dn = Tk.Button(root, text = 'De-noise', command = DeNoise, font=("Copperplate"), fg = 'white',highlightbackground='#0059b3')

# Define slider
HOP_SCALE = Tk.Scale(root,
  length = 200, orient = Tk.HORIZONTAL, from_ = 0.005, to = 3, resolution = 0.005,
  command = hopUpdate,
  label = '     Hop Ratio:',
  font=("Copperplate", 12),
  variable = HOP)
HOP.set(hPr)

VOLUME_SCALE = Tk.Scale(root,
  length = 200, orient = Tk.HORIZONTAL, from_ = 0, to = 10, resolution = 0.01,
  command = VolUpdate,
  label = 'Volume:',
  font=("Copperplate", 12),
  variable = VOL)
VOL.set(gain)

x_pl = 900
y_pl = 120
# Place widgets
B_n.place(height=30, width=150, x=x_pl+25, y= y_pl)
B_r.place(height=30, width=150, x=x_pl+25, y= y_pl+40)
B_w.place(height=30, width=150, x=x_pl+25, y= y_pl+80)
B_t.place(height=30, width=150, x=x_pl+25, y= y_pl+120)
HOP_SCALE.place(x=x_pl, y= y_pl+150)
VOLUME_SCALE.place(x=290, y=430)
B_dn.place(height=30, width=150, x=x_pl+25, y= y_pl+205)
B_quit.place(x=x_pl+70, y= y_pl + 245)
L1.place(x=290, y=15)

#Converting pyplot into widget
canvas = FigureCanvasTkAgg(fig, master = root)
canvas.draw()

W1 = canvas.get_tk_widget()
W1.place(height=350, width=900, x=0, y=60)

#----------------------------------------------------------------------------------------------------
#----------------------------------------READING INPUT-----------------------------------------------
#----------------------------------------------------------------------------------------------------

# Get block of samples from wave file

#----------------------------------------------------------------------------------------------------
#----------------------------------------BASIC ALGORITHM---------------------------------------------
#----------------------------------------------------------------------------------------------------

for i in range(0, num_blocks):
    if QUIT == False:
        root.update()
        #fig.canvas.draw()
        input_bytes = stream.read(BLOCKLEN, exception_on_overflow = False)
        signal_block = struct.unpack('h' * BLOCKLEN, input_bytes)
        
        if a == True:
#----------------------------------------------------------------------------------------------------
#----------------------------------------Whisperization----------------------------------------------
#----------------------------------------------------------------------------------------------------
            # Convert binary data to number sequence (tuple)
            f, t, Zxx = signal.stft(signal_block, RATE, window = 'hann', nperseg=256, noverlap = None)
            lx, ly = Zxx.shape
            phi = numGen(lx, ly)
            whis_abs = np.ndarray(shape=(zxx_len,ly), dtype=complex)
            robot_abs = np.absolute(Zxx)
            whis_abs = robot_abs*phi
            to_be_inverse = whis_abs
            tnew, inver_x = signal.istft(to_be_inverse, RATE, window = 'hann', nperseg=256, noverlap = None)
            inver_x = signal.resample(inver_x, BLOCKLEN)
        elif b == True:
#----------------------------------------------------------------------------------------------------
#----------------------------------------Robotization------------------------------------------------
#----------------------------------------------------------------------------------------------------
            # Convert binary data to number sequence (tuple)
            
            f, t, Zxx = signal.stft(signal_block, 512, window = 'hann', nperseg=stft_N, noverlap = None)
            robot_abs = np.absolute(Zxx)
            to_be_inverse = robot_abs
            tnew, inver_x = signal.istft(to_be_inverse, 512, window = 'hann', nperseg=stft_N, noverlap = None)
            inver_x = signal.resample(inver_x, BLOCKLEN)
            
        elif c == True:
#----------------------------------------------------------------------------------------------------
#--------------------------------------------Normal--------------------------------------------------
#----------------------------------------------------------------------------------------------------
            # Convert binary data to number sequence (tuple)
            f, t, Zxx = signal.stft(signal_block, RATE, window = 'hann', nperseg=stft_N, noverlap = None)
            to_be_inverse = Zxx
            tnew, inver_x = signal.istft(to_be_inverse, RATE, window = 'hann', nperseg=stft_N, noverlap = None)
        elif d == True:
#----------------------------------------------------------------------------------------------------
#-------------------------------------------Pitch-Shift----------------------------------------------
#----------------------------------------------------------------------------------------------------
            # Convert binary data to number sequence (tuple)
            Ni123 = 512
            f, t, Zxx = signal.stft(signal_block, RATE, window = 'hann', nperseg= Ni123, noverlap = Ni123 - ALen)
            Z_mag = np.absolute(Zxx)
            lx, ly = Zxx.shape
            temp = np.ndarray((lx,ly), dtype=complex)
            if flag == 0:
                y_pang = np.ndarray((lx,ly), dtype=complex)
                y_ang = np.ndarray((lx,ly), dtype=complex)
                uwd = np.ndarray((lx,ly), dtype=complex)
                for j in range(0,ly):
                    for i in range(0,lx):
                        uwd[i,j] = (i*2*math.pi*ALen)/Ni123    #frequency * Analysis Length = 2*pi*i* Analysis Length/FFT length
                flag = 1
            y_pang = y_ang      #saving previous angle
            y_ang = np.angle(Zxx)   #getting new data
            yw = (y_ang - y_pang) - uwd;    #(new angle - previous angle) - frequency * Analysis Length
            
            yw = uwd - (np.round_(yw/(2*math.pi))*2*math.pi)
            yw = (yw * hPr)

            if loop1 == True:
                ys = y_ang
                loop1 = False
            else:
                ys = (ys+yw)

            for j in range(ly):
                for i in range(lx):
                   temp[i,j] = Z_mag[i,j]*CmplxAng(np.absolute(ys[i,j]))
            tnew, inver_x = signal.istft(temp, RATE, window = 'hann', nperseg=Ni123, noverlap = Ni123 - SLen)
            inver_x = signal.resample(inver_x, BLOCKLEN)
        elif e == True:
#----------------------------------------------------------------------------------------------------
#----------------------------------------De-Noising--------------------------------------------------
#----------------------------------------------------------------------------------------------------
            # Convert binary data to number sequence (tuple)
            f, t, Zxx = signal.stft(signal_block, RATE, window = 'hann', nperseg=stft_N, noverlap = None)
            Z_mag = np.absolute(Zxx)
            lx, ly = Zxx.shape
            Thres = 0.5
            Z_n_mag = np.ndarray((lx,ly), dtype=complex)
            for j in range(ly):
                    for i in range(lx):
                        if Z_mag[i,j] <= Thres:
                            Zxx[i,j] = 0

            to_be_inverse = Zxx
            tnew, inver_x = signal.istft(to_be_inverse, RATE, window = 'hann', nperseg=stft_N, noverlap = None)
    
        inver_x_int =  np.int_(inver_x)
        g2.set_ydata((inver_x_int*gain)-50)
    
        output_block = np.clip(np.int_(inver_x_int*gain), -MAXVALUE, MAXVALUE)
        output_bytes = struct.pack('h' * BLOCKLEN, *output_block)
    
        # Write binary data to audio output stream
        stream.write(output_bytes, BLOCKLEN)

#----------------------------------------------------------------------------------------------------
#----------------------------------------CLOSING THE CODE--------------------------------------------
#----------------------------------------------------------------------------------------------------
stream.stop_stream()
stream.close()
p.terminate()

pyplot.ioff()           # Turn off interactive mode
pyplot.close()


print('* Finished')

"""coil scan analysis"""

from matplotlib import pyplot as plt
import csv
import numpy as np

plot_type = '2D' # '1D' or '2D'
datafile = "./20230417_164351_coil_scan_zmax_250mV_xmax_ymax_300mV.csv"
has_header = True
with open(datafile, 'r', newline='') as f:
    reader = csv.reader(f)
    if has_header:
        header = reader.__next__()
    data = [row for row in reader]#[:5]
    f.close()

data = np.array(data,float).transpose()
# print(data)
rows,cols = data.shape

if plot_type == '1D':
    # plot all the data in one big plot
    # xpts = range(cols)
    # fig,axes = plt.subplots(nrows=rows,ncols=1,figsize=(20,6))
    # for ax,dat,label in zip(axes,data,header):
    #     ax.plot(xpts,dat)
    #     ax.set_ylabel(label)
    #     ax.set_xlabel("Measurement index")
    # plt.show()

    # split the data into smaller chunks and save plots for each chunk
    chunks = 20
    x_interval = int(cols/chunks) # number of measurements to plot
    xpts = range(x_interval)
    for i in range(chunks):
        xmin = i*x_interval
        xmax = xmin + x_interval
        if xmax > cols-1:
            xmax = -1
        fig,axes = plt.subplots(nrows=rows,ncols=1,figsize=(20,6),dpi=200)
        for ax,dat,label in zip(axes,data,header):
            ax.plot(xpts,dat[xmin:xmax])
            ax.set_ylabel(label)
            ax.set_xlabel("Measurement index")
        plt.savefig(datafile[:-4]+'_'+str(i)+'.png',bbox_inches='tight')
        plt.show()

elif plot_type == '2D':
    counts = data[0]
    zdata,xdata,ydata = data[2:]
    zsteps = len(list(set(zdata))) # figure out how many zsteps there were. this was the outer loop
    xsteps = len(list(set(xdata)))
    ysteps = len(list(set(ydata)))
    sigma = np.std(counts)
    for k in range(zsteps):
        fig, ax = plt.subplots()
        imdata = [[counts[i*xsteps+j+k*xsteps*ysteps] for i in range(xsteps)] for j in range(ysteps)]

        maxcounts = np.amax(imdata)

        cax = ax.imshow(imdata, extent=[min(xdata),max(xdata),max(ydata),min(ydata)])
        ax.set_aspect((max(xdata)-min(xdata))/(max(ydata)-min(ydata)))
        title = "Z volts = {:.3f}".format(zdata[k*xsteps*ysteps])+" V (step "+str(k+1)+'/'+str(zsteps)+')' + \
            '\n Max counts = {}, {:.2f}sigma'.format(maxcounts, maxcounts/sigma)
        ax.set_title(title)
        ax.set_xlabel("X volts")
        ax.set_ylabel("Y volts")

        # test to make sure I understand the axes orientation
        # imdata = np.array([[i*(j+10) for i in range(xsteps)] for j in range(ysteps)])
        # cax = ax.imshow(imdata, extent=[0,xsteps,ysteps+10,10])

        fig.colorbar(cax)
        plt.savefig(datafile[:-4] + '_2D_' + str(k) + '.png', bbox_inches='tight')
        # plt.show()
        plt.close()
else:
    print("sorry,", plot_type," is not a valid plot type. try '1D' or '2D'.")
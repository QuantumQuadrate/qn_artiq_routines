"""coil scan analysis"""

from matplotlib import pyplot as plt
import csv
import numpy as np

plot_type = '3D' # '1D' or '2D' if the scan was 3 dimensional. use '3D' for 4-d data
datafile = "C:\\Networking Experiment\\artiq codes\\artiq-master\\results\\20230504_093318_coil_scan.csv"
has_header = True

if __name__ == '__main__':

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

    elif plot_type == '3D':
        counts = data[0]
        pd_data = data[-1] # the cooling laser power monitor
        coil_data = data[1:-1][:12800] # just as a check. should make two figures

        # sort data so that the zbottom loop which was the inner-most loop is now effectively the second loop
        coil_data = sorted(coil_data,key=lambda x: x[3])
        coil_data = sorted(coil_data,key=lambda x: x[0])

        zbottom_data, ztop_data, xdata, ydata = coil_data

        zbottom_steps = len(list(set(zbottom_data)))
        ztop_steps = len(list(set(zbottom_data)))
        xsteps = len(list(set(xdata)))
        ysteps = len(list(set(ydata)))
        nrows = 4
        ncols = 4
        assert zbottom_steps % nrows == 0, "zbottom_steps should be divisible by nrows for the grid of plots"
        assert ztop_steps % ncols == 0, "ztop_steps should be divisible by ncols for the grid of plots"

        sigma = np.std(counts)
        for k in range(ztop_steps//ncols):
            for l in range(zbottom_steps//nrows):

                fig, axes = plt.subplots(nrows=nrows,ncols=ncols)

                for col in range(ncols):
                    for row in range(nrows):

                        ax = axes[row,col]

                        ztop_step = col*ztop_steps
                        zbottom_steps = row*zbottom_steps
                        imdata = [[counts[i*xsteps + j + k*xsteps*ysteps*col + l*k*xsteps*ysteps*ztop_steps*col*row]
                                   for i in range(xsteps)] for j in range(ysteps)]

                        maxcounts = np.amax(imdata)

                        im = ax.imshow(imdata, extent=[min(xdata), max(xdata), max(ydata), min(ydata)])
                        ax.set_aspect((max(xdata) - min(xdata)) / (max(ydata) - min(ydata)))

                        # this part depends on the ordering of the loops in the scan.
                        # this assumes z_top is the outer-most loop and z_bottom is the inner most
                        # title = "Z_top = {:.3f}V, Z_bottom = {:.3f}V".format(ztop_data[k * xsteps * ysteps * zbottom_steps],
                        #                                                      zbottom_data[l]) \
                        #         + " (step " + str(k + 1) + '/' + str(
                        #     zsteps) + ')' + \
                        #         '\n Max counts = {}, {:.2f}sigma'.format(maxcounts, maxcounts / sigma)
                        # ax.set_title(title)
                        ax.set_xlabel("X volts")
                        ax.set_ylabel("Y volts")

                        # test to make sure I understand the axes orientation
                        # imdata = np.array([[i*(j+10) for i in range(xsteps)] for j in range(ysteps)])
                        # cax = ax.imshow(imdata, extent=[0,xsteps,ysteps+10,10])

                fig.subplots_adjust(right=0.8)
                cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
                fig.colorbar(im, cax=cbar_ax)

                # fig.colorbar(cax)
                # plt.savefig(datafile[:-4] + '_3D_' + str(k) + '.png', bbox_inches='tight')
                plt.show()
                # plt.close()

    else:
        print("sorry,", plot_type," is not a valid plot type.")
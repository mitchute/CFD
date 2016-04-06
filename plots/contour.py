import matplotlib.pyplot as plt
import numpy as np

def createContourPlot( fileName ):

    # Get method name for later
    methodName = fileName.split('.')[0].title()

    # Get data for contour plot
    u = np.genfromtxt(fileName,delimiter=",")

    # Get x locations and delete first entry
    x = u[0]
    x = np.delete(x,[0])

    # Get y locations and delete first entry
    y = u[:,0]
    y = np.delete(y,[0])

    # Get num points
    n = x.shape[0] - 1

    # Remove x location row
    u = np.delete(u,[0],0)

    # Remove y location column
    u = np.delete(u,np.s_[:1],1)

    # Set plot ticks befor running meshgrid
    horizTicks = np.arange(x[0],x[-1]+(x[1] - x[0]),(x[1] - x[0]))
    verTicks = np.arange(y[0],y[-1]+(y[1] - y[0]),(y[1] - y[0]))

    # Configure x and y for plotting
    x, y = np.meshgrid(x, y, indexing='ij')

    # Plot solution
    fig, ax = plt.subplots()
    CS = ax.contour(x, y, u, 45, linewidths=2.5)

    # Tick marks
    ax.set_xticks(horizTicks, minor=True)
    ax.set_yticks(verTicks, minor=True)

    # Grid lines
    ax.grid(which='minor')

    # Title
    plt.title('%s--%d x %d' %(methodName,n,n))

    # Add color bar
    cbar = plt.colorbar(CS, fraction=0.0234, pad=0.04)
    cbar.ax.set_ylabel('$u(x,y)$')

    # Plot formatting
    plt.ylabel('y')
    plt.xlabel('x')

    # Equal aspect ratio
    ax.set_aspect('equal')

    saveFileName = "%s-%d.pdf" %(methodName,n)

    # Save the figure
    plt.savefig(saveFileName, bbox_inches='tight')

    
createContourPlot("explicit.csv")

createContourPlot("implicit.csv")

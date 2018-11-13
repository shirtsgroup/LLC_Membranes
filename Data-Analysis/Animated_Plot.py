import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import math

# Data to be Plotted:
x = [[ 2.77734167,  2.45430833,  2.45983333,  2.44315   ,  2.43258333,
        2.43935   ,  2.43288333,  2.43551667,  2.43945   ,  2.43758333,
        2.43024167,  2.42995833,  2.43750833,  2.44053333,  2.43649167,
        2.43565833,  2.425     ,  2.43394167,  2.442425  ,  2.45350833,
        2.42295833,  2.43779167,  2.43835833,  2.43804167,  2.43186667,
        2.44633333,  2.460475  ,  2.457225  ,  2.44389167,  2.42928333,
        2.4319    ,  2.42885   ,  2.4236    ,  2.43040833,  2.42461667,
        2.42854167,  2.430375  ,  2.44314167,  2.4147    ,  2.424625  ,
        2.4216    ,  2.4268    ,  2.42696667,  2.44165833,  2.44040833,
        2.4308    ,  2.42159167,  2.42176667,  2.42975833,  2.42325   ,
        2.42888333], [-1.26743333, -0.96945833, -0.96101667, -0.95745833, -0.95981667,
       -0.961325  , -0.96300833, -0.974     , -0.97925   , -0.97165833,
       -0.977025  , -0.98663333, -0.97561667, -0.99006667, -0.98300833,
       -0.985675  , -0.98879167, -0.98578333, -0.98263333, -0.97360833,
       -0.98586667, -0.97611667, -0.98363333, -0.99093333, -0.99451667,
       -0.99165   , -0.984     , -0.98295833, -0.98571667, -0.9875    ,
       -0.98018333, -0.982225  , -0.98546667, -0.978025  , -0.980725  ,
       -0.98013333, -0.990525  , -0.98865833, -0.98338333, -0.97296667,
       -0.98166667, -0.981925  , -0.98623333, -0.98525   , -0.99359167,
       -0.99286667, -0.99981667, -0.995475  , -0.98499167, -0.99465833,
       -1.00083333], [ 4.761975  ,  4.31065   ,  4.28734167,  4.26724167,  4.26654167,
        4.26991667,  4.25951667,  4.27685833,  4.26399167,  4.2674    ,
        4.27599167,  4.27263333,  4.272175  ,  4.2612    ,  4.27159167,
        4.25958333,  4.26380833,  4.26849167,  4.26015833,  4.25121667,
        4.25201667,  4.24521667,  4.25225   ,  4.25349167,  4.25663333,
        4.255375  ,  4.25095833,  4.25875833,  4.25240833,  4.26013333,
        4.26163333,  4.25406667,  4.24805   ,  4.25795   ,  4.237275  ,
        4.25265833,  4.240625  ,  4.23125833,  4.23496667,  4.23379167,
        4.23814167,  4.24488333,  4.23448333,  4.24694167,  4.235175  ,
        4.24739167,  4.25015   ,  4.244125  ,  4.23679167,  4.2608    ,
        4.25786667], [ 0.752075  ,  0.82013333,  0.82984167,  0.8281    ,  0.83824167,
        0.83131667,  0.832775  ,  0.84410833,  0.84776667,  0.84695833,
        0.85026667,  0.85745833,  0.8464    ,  0.84448333,  0.84108333,
        0.84928333,  0.828425  ,  0.82455   ,  0.82008333,  0.81573333,
        0.81858333,  0.808225  ,  0.81479167,  0.812925  ,  0.82355   ,
        0.80606667,  0.81349167,  0.80715833,  0.81575833,  0.816575  ,
        0.817225  ,  0.8073    ,  0.817475  ,  0.814375  ,  0.81771667,
        0.8116    ,  0.81275833,  0.81628333,  0.81291667,  0.8226    ,
        0.81395   ,  0.81715   ,  0.82248333,  0.81028333,  0.81741667,
        0.82420833,  0.82105   ,  0.82125   ,  0.81505833,  0.82505833,
        0.81833333]]
y = [[ 4.73530833,  4.14328333,  4.12074167,  4.1637    ,  4.14685833,
        4.15014167,  4.15585833,  4.16618333,  4.1619    ,  4.17358333,
        4.16513333,  4.16206667,  4.174725  ,  4.18005833,  4.16191667,
        4.15340833,  4.12366667,  4.13110833,  4.14231667,  4.13249167,
        4.13964167,  4.14366667,  4.13516667,  4.12994167,  4.13319167,
        4.136375  ,  4.148     ,  4.14016667,  4.13830833,  4.12778333,
        4.13408333,  4.13095   ,  4.1241    ,  4.12710833,  4.12214167,
        4.112075  ,  4.11488333,  4.103825  ,  4.099425  ,  4.10941667,
        4.091775  ,  4.09963333,  4.09925833,  4.09706667,  4.108175  ,
        4.11010833,  4.101425  ,  4.10680833,  4.12056667,  4.11385833,
        4.10161667], [ 4.82683333,  4.37383333,  4.35183333,  4.34098333,  4.33599167,
        4.33293333,  4.318575  ,  4.31903333,  4.30728333,  4.32816667,
        4.3195    ,  4.323875  ,  4.31280833,  4.29910833,  4.30995833,
        4.30991667,  4.32783333,  4.31930833,  4.314575  ,  4.31443333,
        4.31421667,  4.30705833,  4.323025  ,  4.31043333,  4.31721667,
        4.32361667,  4.31805   ,  4.321     ,  4.31879167,  4.317475  ,
        4.32489167,  4.30985   ,  4.30778333,  4.31985833,  4.29524167,
        4.31746667,  4.30144167,  4.288025  ,  4.29879167,  4.28411667,
        4.28890833,  4.29168333,  4.29691667,  4.29898333,  4.29963333,
        4.30598333,  4.30841667,  4.30020833,  4.301175  ,  4.30791667,
        4.29593333], [ 1.234725  ,  1.411925  ,  1.40283333,  1.39491667,  1.394025  ,
        1.39553333,  1.410525  ,  1.4162    ,  1.39938333,  1.40433333,
        1.4089    ,  1.40115   ,  1.40661667,  1.41063333,  1.39938333,
        1.395525  ,  1.40491667,  1.40679167,  1.39803333,  1.38956667,
        1.38195   ,  1.39035   ,  1.38319167,  1.3775    ,  1.390525  ,
        1.38325   ,  1.3818    ,  1.37751667,  1.39365833,  1.366475  ,
        1.38555833,  1.36688333,  1.37079167,  1.37128333,  1.37003333,
        1.36378333,  1.35443333,  1.36668333,  1.357325  ,  1.367225  ,
        1.360075  ,  1.361825  ,  1.36325833,  1.3544    ,  1.36488333,
        1.37659167,  1.37875833,  1.37648333,  1.36741667,  1.36770833,
        1.377425  ], [ 1.340575  ,  1.49430833,  1.511025  ,  1.50783333,  1.51390833,
        1.51105   ,  1.50930833,  1.49750833,  1.50013333,  1.52205833,
        1.51036667,  1.50355   ,  1.50023333,  1.5012    ,  1.50169167,
        1.51156667,  1.49328333,  1.50100833,  1.50494167,  1.514775  ,
        1.51605   ,  1.509175  ,  1.52274167,  1.53078333,  1.525675  ,
        1.53025   ,  1.52065   ,  1.52559167,  1.51810833,  1.5239    ,
        1.51666667,  1.5158    ,  1.49935833,  1.51928333,  1.51365833,
        1.50794167,  1.51298333,  1.50260833,  1.50410833,  1.50980833,
        1.50095   ,  1.50878333,  1.50489167,  1.502325  ,  1.50314167,
        1.50605833,  1.49634167,  1.50458333,  1.51373333,  1.52735   ,
        1.51065833]]

pts, traj_pts = np.shape(x)

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(-2, 6), ylim=(-2, 6))
line, = ax.plot([], [], lw=2)  # Can include data here

# Create an initialization function in order to hide any plot elements that are not wanted to be shown in every frame
def init():
    line.set_data([], [])
    return line,

# Animate function. This is called sequentially. The purpose is to update the plot elements for each frame
def animate(i):
    x_traj = [x[0][i], x[1][i], x[3][i], x[2][i], x[0][i]]  # [Pore1, Pore2, Pore4, Pore3]
    y_traj = [y[0][i], y[1][i], y[3][i], y[2][i], y[0][i]]
    dist = []
    dist.append(math.sqrt((x_traj[0] - x_traj[1])**2 + (y_traj[0] - y_traj[1])**2))  # Pore1 to Pore2
    dist.append(math.sqrt((x_traj[1] - x_traj[2])**2 + (y_traj[1] - y_traj[2])**2))  # Pore2 to Pore4
    dist.append(math.sqrt((x_traj[3] - x_traj[2])**2 + (y_traj[3] - y_traj[2])**2))  # Pore4 to Pore3
    dist.append(math.sqrt((x_traj[0] - x_traj[3])**2 + (y_traj[0] - y_traj[3])**2))  # Pore1 to Pore3
    line.set_data(x_traj, y_traj)
    annotate1 = ax.annotate(('%1.3f A' %dist[0]), xy=((x_traj[0] + x_traj[1])/2, (y_traj[0] + y_traj[1])/2 + .25))
    annotate2 = ax.annotate(('%1.3f A' %dist[1]), xy=((x_traj[2] + x_traj[1])/2 - 1.25, (y_traj[2] + y_traj[1])/2))
    annotate3 = ax.annotate(('%1.3f A' %dist[2]), xy=((x_traj[3] + x_traj[2])/2 - .25, (y_traj[3] + y_traj[2])/2 - .5))
    annotate4 = ax.annotate(('%1.3f A' %dist[3]), xy=((x_traj[0] + x_traj[2])/2 + 2.25, (y_traj[0] + y_traj[2])/2))
    return line, annotate1, annotate2, annotate3, annotate4


# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=traj_pts, interval=2000, blit=True)

plt.title('Pore Center Trajectory')
plt.ylabel('Y Position of Pore Center')
plt.xlabel('X Position of each Pore Center')
plt.show()

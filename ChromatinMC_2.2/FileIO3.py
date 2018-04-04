import matplotlib as mpl

try:
    mpl.use(u'TkAgg')
    mpl.interactive(False)
except:
    pass

import time
from datetime import datetime
import os as os
import numpy as np
import pandas as pd
import sys
from collections import OrderedDict
from etaprogress.progress import ProgressBar
from lmfit import Parameters
from datetime import datetime
import getpass
import distutils.dir_util
import glob
import cv2
from matplotlib import pyplot as plt
from helixmc.pose import HelixPose
# ChromatinMC modules:
import NucleosomeMC3 as nMC
import POVutils as pov

default_folder = 'E:\\users\\'
kT = 41


def _merge(df1, df2):
    r = list(df2.index.values)[0]
    if r in list(df1.index.values):  # existing row
        df3 = df1
        for c in list(df2):
            df3.at[r, c] = df2.at[r, c]
    else:  # new row
        df3 = df1.append(df2)
    return df3


def report_progress(value, title='', init=False):
    global bar, start_time
    if init:
        start_time = time.time()
        print(datetime.now().strftime('>>> Start [{0}] @ %Y-%m-%d %H:%M:%S'.format(title)))
        bar = ProgressBar(value, max_width=100)
    else:
        bar.numerator = value
    elapsed_time = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))

    sys.stdout.write('\r>>> {0} = {1} {2}'.format(elapsed_time, bar, title))
    sys.stdout.flush()
    if bar.done:
        print ('\n')
    return


def get_filename(root=None, ext='dat', incr=False, sub=False, folder=False):
    global working_directory, root_name, file_number, sub_number
    if 'root_name' not in globals():
        if root is None:
            root_name = 'data'
    if root is not None:
        root_name = root
    if 'working_directory' not in globals():
        user = getpass.getuser()
        working_directory = default_folder + '{0:s}\\data\\{1:s}\\'.format(user, datetime.now().strftime('%Y%m%d'))
    if 'file_number' not in globals():
        try:
            filename = working_directory + root_name + '_' + '*.*'
            file_number = int(sorted(glob.glob(filename))[-1].split('_')[-1].split('.')[-2])
        except:
            file_number = 0
    if 'sub_number' not in globals():
        sub_number = 0
    if incr:
        if sub:
            sub_number += 1
        else:
            file_number += 1
            sub_number = 0
    if sub:
        filename = working_directory \
                   + root_name + '_' + str(file_number).zfill(3) + '\\' \
                   + root_name + '_' + str(file_number).zfill(3) + '_' + str(sub_number).zfill(4) + '.' + ext
    else:
        filename = working_directory \
                   + root_name + '_' + str(file_number).zfill(3) + '.' + ext

    distutils.dir_util.mkpath(os.path.dirname(os.path.abspath(filename)))

    if folder:
        filename = os.path.dirname(filename)
    return filename


def change_extension(filename, extension):
    base, _ = os.path.splitext(filename)
    if extension.count('.') == 0:
        extension = '.' + extension
    return base + extension


def read_dataset_xlsx(filename, dataset, pars=None):
    filename = change_extension(filename, 'xlsx')
    sets, files, _ = contents_xlsx(filename)

    row = str(dataset) + ' > ' + files[sets == dataset]

    if pars is None:
        df = pd.read_excel(filename, 'Value')
        pars = Parameters()
        for par in list(df):
            pars.add(str(par), value=df.at[row, par])
        return (pars)

    for p in pars:
        df = pd.read_excel(filename, 'Value')
        try:
            pars[p].value = df.at[row, p]
        except Exception, e:
            pass
        df = pd.read_excel(filename, 'Min')
        try:
            pars[p].min = df.at[row, p]
        except:
            pass
        df = pd.read_excel(filename, 'Max')
        try:
            pars[p].max = df.at[row, p]
        except:
            pass
        df = pd.read_excel(filename, 'StErr')
        try:
            if df.at[row, p] == 0:
                pars[p].vary = False
            else:
                pars[p].vary = True
        except:
            pass
    return pars, files[sets == dataset]


def read_param_xlsx(filename, param):
    filename = change_extension(filename, 'xlsx')
    df = pd.read_excel(filename, 'Value')
    try:
        data= np.asarray(df[param])
        return data
    except:
        print ('Parameter [{:}] not found'.format(param))
    return


def write_xlsx(filename, dataset, pars, report_file=None):
    if not (report_file):
        report_file = filename
    report_file = change_extension(report_file, 'xlsx')
    writer = pd.ExcelWriter(report_file.format('openpyxl'), engine='openpyxl')
    df_index = dataset + ' > ' + filename

    sheets = ['Value', 'StErr', 'Min', 'Max']
    for sheet in sheets:
        p = OrderedDict()
        if sheet == 'StErr':
            for k, v in pars.items():
                try:
                    p[k] = v.ster
                except:
                    if v.vary:
                        p[k] = np.inf
                    else:
                        p[k] = 0
        elif sheet == 'Min':
            for k, v in pars.items():
                p[k] = v.min
        elif sheet == 'Max':
            for k, v in pars.items():
                p[k] = v.max
        elif sheet == 'Value':
            for k, v in pars.items():
                p[k] = v.value

        df = pd.DataFrame(p, index=[df_index])
        try:
            old_dfs = pd.read_excel(report_file, sheet)
            df = _merge(old_dfs, df)
        except:
            pass
        df.to_excel(writer, sheet_name=sheet)

    writer.close()
    return df_index


def contents_xlsx(filename):
    filename = change_extension(filename, 'xlsx')
    df = pd.read_excel(filename, 'Value')
    datasets = []
    files = []
    for txt in list(df.index.values):
        datasets.append(int(txt.split(' > ')[0]))
        files.append(str(txt.split(' > ')[1]))
        par_names = []
        for item in list(df):
            par_names.append(str(item))
    return datasets, files, par_names


def create_mp4_images(image_files, fps=5.0, delete=True, filename=None, reverse=True):
    image_files = sorted(image_files)
    if filename == None:
        filename = image_files[0]
    filename = change_extension(filename, 'mp4')

    frame = cv2.imread(image_files[0])
    height, width, channels = frame.shape

    # fourcc = cv2.VideoWriter_fourcc(*'mp4v') # Be sure to use lower case
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')  # Be sure to use lower case
    out = cv2.VideoWriter(filename, fourcc, fps, (width, height))

    for image in image_files:
        frame = cv2.imread(image)
        out.write(frame)  # Write out frame to video
    if reverse:
        for image in reversed(image_files):
            frame = cv2.imread(image)
            out.write(frame)  # Write out frame to video
    if delete:
        for image in image_files:
            os.remove(image)
    out.release()
    print('Movie saved as:', filename)
    return filename


def save_plot(data, ax_labels=None, grid=None, xrange=None, yrange=None, save=True, filename=None, transpose=False):
    #   data = [x, y, line, selected]
    if filename is None:
        filename = get_filename(incr=True)
    plt.close()
    plt.figure(figsize=(4, 3))
    if not transpose:
        xsel = data[0][data[3] != 0]
        ysel = data[1][data[3] != 0]
        xunsel = data[0][data[3] == 0]
        yunsel = data[1][data[3] == 0]
        xline = data[0]
        yline = data[2]
    else:
        xsel = data[1][data[3] != 0]
        ysel = data[0][data[3] != 0]
        xunsel = data[1][data[3] == 0]
        yunsel = data[0][data[3] == 0]
        xline = data[2]
        yline = data[0]
    plt.scatter(xunsel, yunsel, s=10, facecolors='none', edgecolors='grey')
    plt.scatter(xsel, ysel, s=10, facecolors='none', edgecolors='b')
    plt.plot(xline, yline, color='k')

    if grid is not None:
        for g in grid:
            if not transpose:
                plt.plot(xline, g, color='k', linestyle=':', linewidth=0.8)
            else:
                plt.plot(g, yline, color='k', linestyle=':', linewidth=0.8)

    if ax_labels is not None:
        plt.xlabel(ax_labels[0])
        plt.ylabel(ax_labels[1])
    if xrange is not None:
        plt.xlim(xrange)
    if yrange is not None:
        plt.xlim(yrange)
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    plt.title(filename.split('\\')[-2] + '\\' + filename.split('\\')[-1].split('.')[0], loc='left',
              fontdict={'fontsize': 10})
    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.draw()
    plt.pause(0.05)
    if save:
        filename = get_filename(ext='png')
        plt.savefig(filename, dpi=600, format='png')
    return filename


def pov_dna(filename, coords, range_A=[1000, 1000], offset_A=[0, 0, 0], width_pix=500, colors=None, radius=None,
            show=False):
    if radius is None:
        radius = np.append([10], (np.ones(8) * 4))
        radius = np.append(radius, 3)
        radius = np.append(radius, 2)
    if colors is None:
        colors = 'kbbggryrycy'

    filename = change_extension(filename, 'pov')

    aspect_ratio = range_A[1] / float(range_A[0])
    pov_image = pov.init(plt_width=range_A[0], aspect_ratio=aspect_ratio)
    i = 0
    j = 0
    offset = np.asarray(offset_A) - np.asarray((0, 0, range_A[1] / 2.0))
    for coord in coords:
        if (i > len(colors) - 1):
            i = 0
        if (j > len(radius) - 1):
            j = 0
        for sphere in coord:
            pov_image = pov.add_sphere(pov_image, sphere + offset, color=colors[i], radius=radius[j])
        i += 1
        j += 1
    pov.save(pov_image, filename=filename)
    pov.render(filename, height=width_pix * aspect_ratio, width=width_pix)
    if show:
        pov.show(filename)
    filename = change_extension(filename, 'png')
    # print('File created: {:s}'.format(filename))
    return filename


def plot_dna(pose, tf=None, color='blue', update=False, title='', range_nm=100, save=False, wait=0):
    global ax, fig, scale

    if tf is None:
        coords = pose.coords / 10.0
    else:
        coords = nMC.apply_transformation(pose.coords, tf) / 10.0

    if update:
        ax.clear()
        plt.title(title, loc='left')
    else:
        plt.close()
        fig = plt.figure(figsize=(5, 5))
        ax = fig.add_subplot(111, projection='3d')
        scale = range_nm / 2
        plt.title(title, loc='left')
        plt.tight_layout(pad=0.3, w_pad=0.3, h_pad=0.3)
    pointsize = 0.1 * np.sqrt(scale)

    ax.set_xlabel('x (nm)')
    ax.set_ylabel('y (nm)')
    ax.set_zlabel('z (nm)')

    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], s=pointsize, c=color, alpha=0.25)
    ax.scatter(coords[0, 0], coords[0, 1], coords[0, 2], s=pointsize, c='red', alpha=0.55)
    ax.scatter([0], [0], [0], s=pointsize, c='k', alpha=0.55)

    plt.xlim((-scale, scale))
    plt.ylim((-scale, scale))
    ax.set_zlim(0, 2 * scale)

    plt.draw()
    plt.pause(0.05)
    if save:
        filename = get_filename(incr=True, sub=True, ext='jpg')
        plt.savefig(filename, dpi=600, format='jpg')
    time.sleep(wait)
    return


def create_mp4_pov(directory, origin_frame=0, fps=5, reverse=True):
    filenames = sorted(glob.glob(directory + '\\*.npz'))
    r = np.append([10], (np.ones(8) * 4))
    r = np.append(r, 3)
    r = np.append(r, 2)
    c = 'kbbggryrycy'

    image_files = []
    j = 0
    report_progress(len(filenames), title='create_pov_mp4', init=True)
    for file in filenames:
        j += 1
        report_progress(j, title='create_pov_mp4')

        dna = HelixPose.from_file(file)
        origin_of = nMC.get_of(dna, origin_frame)
        tf = nMC.get_transformation(origin_of)
        coords = nMC.apply_transformation(dna.coords, tf)
        image_files.append(pov_dna(file, [coords], range_A=[1000, 1500], offset_A=[0, 0, 150], show=False))
    create_mp4_images(image_files, fps=fps, delete=True, filename=directory + '.mp4', reverse=reverse)
    return



def create_pov_npz(filename):
    filename = change_extension(filename, 'npz')
    # directory = 'E:\\Users\\Noort\\data\\20180320\\data_001\\'
    # filename = 'E:\\Users\\Noort\\data\\20180320\\data_001\\data_001_0001.npz'
    # print(filename)
    dna = HelixPose.from_file(filename)
    pov_dna(filename, [dna.coords], range_A=[1500, 2500], offset_A=[0, 0, 150], show=True)
    return

    # filename = 'E:\\users\\noort\\data\\20180320\\12nucs_016.dat'
    # plot_step_params(filename, -1, save=True)


def main():
    filename = 'E:\\users\\noort\\data\\20180321\\12nucs_003.dat'
    analyze_step_params(filename)
    return


if __name__ == "__main__":
    # execute only if run as a script
    main()
    # create_pov_mp4('E:\\Users\\Noort\\data\\20180315\\2nucs_001', origin_frame=0, reverse=False)
    # filename = get_filename(incr=True, base='2nucs')
    # print(filename)

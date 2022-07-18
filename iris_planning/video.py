import typing as typ
import pathlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.cm
import astropy.units as u
import kgpy.observatories.iris


def convert_folder(
        folder_input: pathlib.Path,
        file_output: pathlib.Path,
        frame_interval: u.Quantity = 1 / 30 * u.s,
        num_files: typ.Optional[int] = None,
        channel_sji: int = 1400,
):
    sji_path_array = np.array(sorted(folder_input.glob(f'*SJI_{channel_sji}*.fits*')))[..., np.newaxis]
    sg_archive_sequence = sorted(folder_input.glob('*raster*'))

    if num_files is not None:
        sji_path_array = sji_path_array[:num_files]
        sg_archive_sequence = sg_archive_sequence[:num_files]

    convert_files(
        sji_path_array=sji_path_array,
        sg_archive_sequence=sg_archive_sequence,
        file_output=file_output,
        frame_interval=frame_interval
    )


def convert_files(
    sji_path_array: np.ndarray,
    sg_archive_sequence: typ.Sequence[pathlib.Path],
    file_output: pathlib.Path,
    frame_interval: u.Quantity = 1 / 30 * u.s,
):

    cube_sequence = []
    for archive in sg_archive_sequence:
        if archive.suffix == '.fits':
            cube = kgpy.observatories.iris.spectrograph.Cube.from_path_sequence([archive])
        elif archive.suffix == '.gz':
            cube = kgpy.observatories.iris.spectrograph.Cube.from_archive(archive)
        else:
            continue
        cube_sequence.append(cube)

    sg_cubes = kgpy.observatories.iris.spectrograph.CubeList(
        cube_sequence
    )

    # sg_cubes = kgpy.observatories.iris.spectrograph.CubeList(
    #     [kgpy.observatories.iris.spectrograph.Cube.from_archive(archive) for archive in sg_archive_sequence]
    # )
    shape_w_min = np.array([c.shape[~0] for c in sg_cubes]).min()
    for c in sg_cubes:
        c.intensity = c.intensity[..., :shape_w_min]

    obs = kgpy.observatories.iris.Obs(
        sji_images=kgpy.observatories.iris.sji.ImageList(
            [kgpy.observatories.iris.sji.Image.from_path_sequence(paths) for paths in sji_path_array]
        ).to_image(),
        sg_cubes=sg_cubes.to_cube(),
    )

    fig, ax = plt.subplots(figsize=(16 / 2, 9 / 2), constrained_layout=True)
    fig.set_facecolor('black')
    ax.set_facecolor('black')
    ani = obs.animate_channel_color(
        ax=ax,
        frame_interval=frame_interval,
        thresh_max=99 * u.percent,
        thresh_min=1 * u.percent,
    )

    ani.save(
        filename=file_output,
        dpi=200,
        progress_callback=lambda i, n: print(f'Saving frame {i} of {n}')
    )


def convert_folder_sg_colors(
        folder_input: pathlib.Path,
        file_output: pathlib.Path,
        frame_interval: u.Quantity = 1 / 30 * u.s,
        max_doppler_shift: u.Quantity = 50 * u.km / u.s,
):
    archive_sequence = sorted(folder_input.glob('*raster*'))

    cube_sequence = []
    for archive in archive_sequence:
        if archive.suffix == '.fits':
            cube = kgpy.observatories.iris.spectrograph.Cube.from_path_sequence([archive])
        elif archive.suffix == '.gz':
            cube = kgpy.observatories.iris.spectrograph.Cube.from_archive(archive)
        else:
            continue
        cube_sequence.append(cube)

    cube = kgpy.observatories.iris.spectrograph.CubeList(
        cube_sequence
        # [kgpy.observatories.iris.spectrograph.Cube.from_archive(archive) for archive in archive_sequence]
    )
    shape_w_min = np.array([c.shape[~0] for c in cube]).min()
    for c in cube:
        c.intensity = c.intensity[..., :shape_w_min]
        # c.intensity_uncertainty = c.intensity_uncertainty[..., :shape_w_min]

    cube = cube.to_cube()

    # cube.intensity = np.nan_to_num(cube.intensity)
    # cube.intensity = cube.intensity_despiked

    fig, ax = plt.subplots(figsize=(8, 8), constrained_layout=True)
    # fig.set_facecolor('black')
    ax.set_facecolor('black')
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    # ax.set_axis_off()

    # max_doppler_shift = 50 * u.km / u.s
    mappable = cube.colormap_spectrum
    mappable = matplotlib.cm.ScalarMappable(
        cmap=mappable.cmap,
        norm=matplotlib.colors.Normalize(vmin=-max_doppler_shift.value, vmax=max_doppler_shift.value)
    )
    colorbar = plt.colorbar(
        mappable=mappable,
        ax=ax,
    )
    colorbar.ax.set_ylabel(max_doppler_shift.unit)
    ani = cube.animate_colors(
        ax=ax,
        frame_interval=frame_interval,
        thresh_max=99 * u.percent,
        thresh_min=1 * u.percent,
        max_doppler_shift=max_doppler_shift,
    )

    # plt.show()

    ani.save(
        filename=file_output,
        dpi=200,
        progress_callback=lambda i, n: print(f'Saving frame {i} of {n}'),
    )

    # import matplotlib as mpl
    # print('bitrate', mpl.rcParams['animation.bitrate'])
    # mpl.rcParams['animation.bitrate'] = 1000000000
    # print('bitrate', mpl.rcParams['animation.bitrate'])
    #
    # with open(file_output.parent / f'{file_output.name}.html', "w") as f:
    #     print(ani.to_html5_video(embed_limit=100), file=f)
    #     # print(ani.to_jshtml(), file=f)



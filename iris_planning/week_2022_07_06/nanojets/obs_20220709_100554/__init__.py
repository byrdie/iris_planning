import pathlib
import astropy.units as u
import iris_planning.video

__all__ = ['to_video']


path_base = pathlib.Path(__file__).parent


def to_video():
    iris_planning.video.convert_folder(
        folder_input=path_base / 'data',
        file_output=path_base / f'obs_20220709_100554.mp4',
        frame_interval=1/24 * u.s,
        channel_sji=1330,
        # num_files=3,
    )


def to_video_sg_colors():
    iris_planning.video.convert_folder_sg_colors(
        folder_input=path_base / 'data',
        file_output=path_base / f'obs_20220709_100554_sg_colors.mp4',
        frame_interval=1/4 * u.s,
    )

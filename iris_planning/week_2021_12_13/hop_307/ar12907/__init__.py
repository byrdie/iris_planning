import pathlib
import astropy.units as u
import iris_planning.video

__all__ = ['to_video']


path_base = pathlib.Path(__file__).parent


def to_video():
    iris_planning.video.convert_folder(
        folder_input=path_base / 'data',
        file_output=path_base / f'ar12907.mp4',
        frame_interval=1/60 * u.s,
        # num_files=3,
    )


def to_video_sg_colors():
    iris_planning.video.convert_folder_sg_colors(
        folder_input=path_base / 'data',
        file_output=path_base / f'ar12907_sg_colors.mp4',
        frame_interval=1/4 * u.s,
    )
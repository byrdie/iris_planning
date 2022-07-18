from . import to_video, to_video_sg_colors


def test_to_video(capsys):
    with capsys.disabled():
        to_video()


def test_to_video_sg_colors(capsys):
    with capsys.disabled():
        to_video_sg_colors()

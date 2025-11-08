# RnaThermofinder/utils/file_converters.py

"""File format conversion utilities"""

from pathlib import Path
from typing import Optional, List
import cairosvg


def svg_to_png(svg_path: Path, png_path: Optional[Path] = None,
               scale: float = 2.0) -> Optional[Path]:
    """
    Convert SVG file to PNG

    Args:
        svg_path: Path to SVG file
        png_path: Path to save PNG (if None, replaces .svg with .png)
        scale: Scale factor for resolution (2.0 = double resolution)

    Returns:
        Path to PNG file if successful, None otherwise

    Example:
       # >>> svg_to_png(Path("structure.svg"))
        Path("structure.png")
    """
    if png_path is None:
        png_path = svg_path.with_suffix('.png')

    try:
        cairosvg.svg2png(
            url=str(svg_path),
            write_to=str(png_path),
            scale=scale
        )
        return png_path
    except Exception as e:
        print(f"Error converting {svg_path.name}: {e}")
        return None


def svg_to_png_batch(svg_files: List[Path], output_dir: Optional[Path] = None,
                     scale: float = 2.0) -> List[Path]:
    """
    Convert multiple SVG files to PNG

    Args:
        svg_files: List of SVG file paths
        output_dir: Directory for PNG files (if None, same as SVG location)
        scale: Scale factor for resolution

    Returns:
        List of successfully created PNG file paths
    """
    png_files = []

    for svg_path in svg_files:
        if output_dir:
            png_path = output_dir / svg_path.with_suffix('.png').name
        else:
            png_path = svg_path.with_suffix('.png')

        result = svg_to_png(svg_path, png_path, scale)
        if result:
            png_files.append(result)

    return png_files


def create_thumbnail(image_path: Path, thumbnail_path: Optional[Path] = None,
                     max_size: tuple = (200, 200)) -> Optional[Path]:
    """
    Create thumbnail from image

    Args:
        image_path: Path to image file
        thumbnail_path: Path to save thumbnail
        max_size: Maximum dimensions (width, height)

    Returns:
        Path to thumbnail if successful
    """
    try:
        from PIL import Image

        if thumbnail_path is None:
            thumbnail_path = image_path.parent / f"{image_path.stem}_thumb{image_path.suffix}"

        img = Image.open(image_path)
        img.thumbnail(max_size, Image.Resampling.LANCZOS)
        img.save(thumbnail_path, optimize=True)
        return thumbnail_path
    except Exception as e:
        print(f"Error creating thumbnail: {e}")
        return None
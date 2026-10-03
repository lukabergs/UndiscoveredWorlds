"""Per-product export geometry; absent metadata preserves historical contracts.

Current runs write maps/map_render_metadata.json (version 2): every PNG map has
the same width. Older runs wrote wind_render_metadata.json (version 1), where
only wind products used that width and scalar maps kept the world grid.
"""
import json

WIND_PRODUCTS = {'lic', 'particles', 'speed'}

def _checked(path, data, version, world_width):
    width = data.get('columns')
    if (data.get('version') != version or type(width) is not int or not 64 <= width <= 2048 or width % 2
            or data.get('rows') != width // 2 or data.get('source_columns') != world_width
            or data.get('source_rows') != world_width // 2):
        raise ValueError(f'Invalid map render geometry: {path}')
    return data

def load_wind_geometry(directory, world_width):
    path = directory / 'maps/map_render_metadata.json'
    if path.exists():
        data = _checked(path, json.loads(path.read_text(encoding='utf-8')), 2, world_width)
        if data.get('products') != 'all':
            raise ValueError(f'Invalid map render geometry: {path}')
        return data
    path = directory / 'maps/wind_render_metadata.json'
    if not path.exists():
        return None
    data = _checked(path, json.loads(path.read_text(encoding='utf-8')), 1, world_width)
    if set(data.get('products', [])) != WIND_PRODUCTS:
        raise ValueError(f'Invalid wind render geometry: {path}')
    return data

def expected_map_size(relative_path, world_width, wind_geometry):
    parts = relative_path.parts
    wind = len(parts) >= 2 and parts[0] == 'wind' and parts[1] in WIND_PRODUCTS
    shared = wind_geometry is not None and wind_geometry.get('products') == 'all'
    width = wind_geometry['columns'] if wind_geometry is not None and (shared or wind) else world_width
    return width, width // 2

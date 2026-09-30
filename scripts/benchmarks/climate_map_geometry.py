"""Per-product export geometry; absent metadata preserves historical contracts."""
import json

WIND_PRODUCTS = {'lic', 'particles', 'speed'}

def load_wind_geometry(directory, world_width):
    path = directory / 'maps/wind_render_metadata.json'
    if not path.exists():
        return None
    data = json.loads(path.read_text(encoding='utf-8'))
    width = data.get('columns')
    if (data.get('version') != 1 or type(width) is not int or not 64 <= width <= 2048 or width % 2
            or data.get('rows') != width // 2 or data.get('source_columns') != world_width
            or data.get('source_rows') != world_width // 2 or set(data.get('products', [])) != WIND_PRODUCTS):
        raise ValueError(f'Invalid wind render geometry: {path}')
    return data

def expected_map_size(relative_path, world_width, wind_geometry):
    parts = relative_path.parts
    wind = len(parts) >= 2 and parts[0] == 'wind' and parts[1] in WIND_PRODUCTS
    width = wind_geometry['columns'] if wind and wind_geometry is not None else world_width
    return width, width // 2

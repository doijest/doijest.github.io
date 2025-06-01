import os
import json

# Paths
input_geojson = "moorwalkers.geojson"
output_dir = "tracks"

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Load GeoJSON
with open(input_geojson, "r", encoding="utf-8") as f:
    data = json.load(f)

features = data.get("features", [])

manifest = []
marker_features = []
for feature in features:
    name = feature.get("properties", {}).get("name")
    if not name:
        continue  # Skip features without a name

    # Sanitize filename
    safe_name = "".join(c for c in name if c.isalnum() or c in (' ', '_', '-')).rstrip()
    filename = f"{safe_name}.geojson"
    filepath = os.path.join(output_dir, filename)

    # Create single-feature GeoJSON
    feature_geojson = {
        "type": "FeatureCollection",
        "features": [feature]
    }

    with open(filepath, "w", encoding="utf-8") as out_f:
        json.dump(feature_geojson, out_f, ensure_ascii=False, indent=2)

    # Add relative path to manifest
    manifest.append(os.path.join(output_dir, filename).replace("\\", "/"))

    # Extract starting location for marker
    geometry = feature.get("geometry", {})
    coords = geometry.get("coordinates", [])
    if geometry.get("type") == "LineString" and coords:
        start_coord = coords[0]
        # Only use lon, lat (ignore elevation/time if present)
        marker_feature = {
            "type": "Feature",
            "geometry": {
                "type": "Point",
                "coordinates": start_coord[:2]
            },
            "properties": feature.get("properties", {})
        }
        marker_features.append(marker_feature)

# Write manifest file
manifest_path = "tracks-manifest.json"
with open(manifest_path, "w", encoding="utf-8") as mf:
    json.dump(manifest, mf, ensure_ascii=False, indent=2)

# Write track_markers.geojson file
track_markers_geojson = {
    "type": "FeatureCollection",
    "features": marker_features
}
with open("track_markers.geojson", "w", encoding="utf-8") as mf:
    json.dump(track_markers_geojson, mf, ensure_ascii=False, indent=2)

print(f"Split {len(features)} features into '{output_dir}' folder, created manifest '{manifest_path}', and created 'track_markers.geojson' with {len(marker_features)} markers.")
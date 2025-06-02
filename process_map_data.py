import os
from datetime import timedelta  # Date and time module
from datetime import datetime
import json
import geojson
import gpxpy
import requests
import time
from geopy import distance  # Geographical distance calculation module
from shapely.geometry import LineString  # Geometric operations library
from OSGridConverter import (
    latlong2grid,
)  # Latitude and longitude to grid reference conversion
from sklearn.cluster import KMeans  # K-means clustering module from scikit-learn
from xml.etree.ElementTree import (
    Element,
    SubElement,
    tostring,
)  # XML data manipulation module
import matplotlib.pyplot as plt  # Visualization library for creating plots and graphs

def find_place_name(starting_latitude, starting_longitude):
    # Construct the Nominatim API URL for reverse geocoding
    nominatim_url = f"https://nominatim.openstreetmap.org/reverse?format=json&lat={starting_latitude}&lon={starting_longitude}"
    response = requests.get(nominatim_url, timeout=10)  # Set a timeout for the request
    if response.status_code == 200:
        data = response.json()
        if data:
            # Extract the nearest place name
            display_name = data.get("display_name", "Unknown")
            # Split the string by commas and join the required parts
            display_name_parts = display_name.split(",")
            if len(display_name_parts) > 5:
                display_name = ",".join(display_name_parts[:-5])
            return display_name
    # If request fails, provide the user with the option to manually input the location name
    print(f"Unable to retrieve data from {nominatim_url}, please set manually.")
    print(f"You can use the following URL to find the display_name: {nominatim_url}")
    print("Please enter the full display_name value:")
    display_name = input()
    # Split the string by commas and join the required parts
    display_name_parts = display_name.split(",")
    if len(display_name_parts) > 5:
        display_name = ",".join(display_name_parts[:-5])
    return display_name

def create_data(main_geojson):
    # Create an empty list to store the track years
    years = []

    # Load existing feature_collection data from previous version of geoJSON file
    # Load the GeoJSON file
    if os.path.exists(main_geojson):
        with open(main_geojson, "r", encoding="utf-8") as f:
            geojson_data = geojson.load(f)
        # Check if it's a FeatureCollection, and if so load it
        if isinstance(geojson_data, geojson.FeatureCollection):
            feature_collection = geojson_data
            print("Existing data loaded")
        else:
            print(
                "The loaded data is not a FeatureCollection. Creating an empty features list."
            )
            features = []
        # Create features list and load years list from the feature_collection
        features = feature_collection["features"]
        for feature in features:
            years.append(feature["properties"]["name"][:4])

    else:
        print("The GeoJSON file does not exist. Creating an empty features list.")
        features = []

    # Create list of existing tracks already loaded as features
    existing_track_list = []
    for feature in features:
        existing_track_list.append(feature["properties"]["name"] + ".gpx")

    # Loop through each file in the input directory
    input_dir = os.path.join(os.getcwd(), "orig_gpx_files")
    file_list = os.listdir(input_dir)
    # Filter out files that exist in existing_track_list
    file_list = [
        filename
        for filename in file_list
        if filename.endswith(".gpx") and filename not in existing_track_list
    ]
    file_list = list(reversed(file_list))
    file_list_length = len(file_list)
    processed_count = 0
    print(f"Processing {file_list_length} new files")

    for filename in file_list:
        # Add the year from the filename to the years list
        years.append(filename[:4])

        # Load the GPX file and parse the data
        with open(
            os.path.join(input_dir, filename), "r", encoding="utf-8"
        ) as file_path:
            gpx_data = gpxpy.parse(file_path)

        # Extract the track segments and points from the GPX file
        track_segments = [s for s in gpx_data.tracks[0].segments]
        track_points = [p for s in track_segments for p in s.points]

        # Calculate cumulative distances for each point, required for elevation profiles
        total_distance = 0.0
        cumulative_distances = [0.0]  # Starting with 0 distance
        for i in range(1, len(track_points)):
            coord1 = (track_points[i - 1].latitude, track_points[i - 1].longitude)
            coord2 = (track_points[i].latitude, track_points[i].longitude)
            dist = distance.distance(coord1, coord2).miles
            total_distance += dist
            cumulative_distances.append(total_distance)

        # Create a GeoJSON feature LineString with track points and cumulative distance for each point
        line_coords = [
            (p.longitude, p.latitude, p.elevation, cumulative_distances[i])
            for i, p in enumerate(track_points)
        ]

        # Calculate the total distance of the track in km
        total_distance_km = 0
        for i in range(len(track_points) - 1):
            total_distance_km += distance.distance(
                (track_points[i].latitude, track_points[i].longitude),
                (track_points[i + 1].latitude, track_points[i + 1].longitude),
            ).km

        # Calculate the total distance of the track in miles
        total_distance_mi = 0
        for i in range(len(track_points) - 1):
            total_distance_mi += distance.distance(
                (track_points[i].latitude, track_points[i].longitude),
                (track_points[i + 1].latitude, track_points[i + 1].longitude),
            ).miles

        # Determine the overall ascent and descent, measured in meters, across the track by utilizing a 7-Point averaged list of elevations
        # Initialize an empty list to store computed values
        sevenpoint = []
        # Calculate the last index of the track_points list
        max_index = len(track_points) - 1  # Subtract 1 to get the last index

        # Iterate through track_points using index and value
        for idx, val in enumerate(track_points):
            # Append elevations at the edges directly to sevenpoint list
            if idx == 0 or idx == max_index:
                sevenpoint.append(track_points[idx].elevation)
            # Compute the average elevation based on the immediate neighbors for second and second-to-last elements
            elif idx == 1 or idx == max_index - 1:
                sum =   track_points[idx].elevation + \
                        track_points[idx-1].elevation + \
                        track_points[idx+1].elevation
                sevenpoint.append(sum / 3)
            # Compute the average elevation based on a window of 5 elements for third and third-to-last elements
            elif idx == 2 or idx == max_index - 2:
                sum =   track_points[idx].elevation + \
                        track_points[idx-1].elevation + \
                        track_points[idx-2].elevation + \
                        track_points[idx+1].elevation + \
                        track_points[idx+2].elevation
                sevenpoint.append(sum / 5)
            # Compute the average elevation based on a window of 7 elements for the rest of the elements
            else:
                sum =   track_points[idx].elevation + \
                        track_points[idx - 1].elevation + \
                        track_points[idx - 2].elevation + \
                        track_points[idx - 3].elevation + \
                        track_points[idx + 1].elevation + \
                        track_points[idx + 2].elevation + \
                        track_points[idx + 3].elevation
                sevenpoint.append(sum / 7)

        # Initialize variables to calculate total ascent and descent
        total_ascent_m = 0
        total_descent_m = 0
        # Iterate through computed elevations in sevenpoint to calculate ascent and descent
        for idx, val in enumerate(sevenpoint):
            # Skip the first index since there's no previous elevation to compare
            if idx != 0:
                elevation_change = sevenpoint[idx] - sevenpoint[idx-1]
                # Calculate total ascent and descent based on elevation change
                if elevation_change > 0:
                    total_ascent_m += elevation_change
                elif elevation_change < 0:
                    total_descent_m += elevation_change

        # Calculate the starting and centre coordinates
        # Extract latitude and longitude from the track points
        coordinates = [(point.latitude, point.longitude) for point in track_points]
        # Create a LineString object from the track points
        line = LineString(coordinates)
        # Calculate the center and starting point of the LineString
        center_point = line.centroid
        starting_point = line.coords[0]
        # Extract latitude and longitude from the center point
        center_latitude = center_point.x
        center_longitude = center_point.y
        # Extract latitude and longitude from the starting point
        starting_latitude = starting_point[0]
        starting_longitude = starting_point[1]

        # Calculate the OS Grid Ref
        # Convert latitude and longitude to grid reference
        gridref = str(latlong2grid(starting_latitude, starting_longitude))
        # print(gridref)

        # Generate the Google Maps link to the starting point
        googleMapsLink = (
            "https://www.google.com/maps?q="
            + str(starting_latitude)
            + ","
            + str(starting_longitude)
        )
        # print(googleMapsLink)

        # Generate the track_download link
        download_link = (
            "https://moorwalkers.github.io/track_downloads/"
            + filename[:-4].replace(" ", "").replace("@", "_")
            + ".gpx"
        )
        # print(download_link)

        # Generate the ind_map links
        ind_map_link = (
            "https://moorwalkers.github.io/ind_maps/"
            + filename[:-4].replace(" ", "").replace("@", "_")
            + ".html"
        )
        # print(ind_map_link)
        ind_map_link_os = (
            "https://moorwalkers.github.io/ind_maps_os/"
            + filename[:-4].replace(" ", "").replace("@", "_")
            + "_os.html"
        )
        # print(ind_map_link_os)

        # Generate the elevation profile image link
        elevation_profile_link = (
            "https://moorwalkers.github.io/elevation_profiles/"
            + filename[:-4].replace(" ", "").replace("@", "_")
            + ".png"
        )
        #print(elevation_profile_link)

        # Calculate the start and end times of the track, and then the duration
        start_time = track_points[0].time
        end_time = track_points[-1].time
        duration = end_time - start_time
        duration = timedelta(seconds=duration.seconds)  # Remove milliseconds

        # Calculate a JavaScript compatible datetime
        if filename[:4] == "2020":
            date_string = "2020-07-01 @ 00-00-00"
        else:
            date_string = filename[:-4]
        # Replace @ with T and - with :
        formatted_date_string = date_string.replace('@', 'T').replace('-', ' ')
        # Adjust the format for parsing the date string
        python_date = datetime.strptime(formatted_date_string, "%Y %m %d T %H %M %S")
        # Convert to ISO 8601 format
        iso_date = python_date.isoformat()
        #print(iso_date)

        # Create a GeoJSON feature LineString from the GPX file
        feature = geojson.Feature(
            geometry=geojson.LineString(line_coords),
            properties={
                "name": filename[:-4],
                "date": iso_date,
                "distance_km": round(total_distance_km, 2),
                "distance_mi": round(total_distance_mi, 2),
                "duration": str(duration),
                "ascent": int(total_ascent_m),
                "descent": int(total_descent_m),
                "center_lat": center_latitude,
                "center_lon": center_longitude,
                "place_name": find_place_name(starting_latitude, starting_longitude),
                "gridref": gridref,
                "googleMapsLink": googleMapsLink,
                "download_link": download_link,
                "ind_map_link": ind_map_link,
                "ind_map_link_os": ind_map_link_os,
                "elevation_profile_link": elevation_profile_link,
            },
        )

        # Add the feature to the features list
        features.append(feature)

        processed_count += 1
        processed_percent = "{:.1f}%".format((processed_count / file_list_length) * 100)
        print(f"{processed_percent} processed")

    # Sort the features list by name
    features = sorted(features, key=lambda x: x["properties"]["name"], reverse=True)

    # Create a GeoJSON feature collection from the features list
    feature_collection = geojson.FeatureCollection(features)

    # Cluster the GeoJSON data into groups of nearby start points
    # Required to decrease the likelihood of the same colour being used for nearby tracks

    # List of colours to be used for the tracks and markers
    colours = [
        "darkred",
        "green",
        "red",
        "blue",
        "gray",
        "purple",
        "black",
        "cadetblue",
        "darkgreen",
        "orange",
        "darkblue",
    ]

    # Calulcate how many clusters should be produced by dividing the total number of tracks by the
    # total number of available colours and adding one, this should ensure no cluster has more
    # tracks than available colours
    num_clusters = len(feature_collection["features"]) // len(colours) + 1
    # print(f"Number of featues = {len(feature_collection['features'])}")
    # print(f"Number of colours = {len(colours)}")
    # print(f"Number of clusters = {num_clusters}")

    # Extract the start point coordinates from GeoJSON file
    coordinates = []
    for feature in feature_collection["features"]:
        # Use the first point in each linestring, just returning the lat and long (not elevation)
        first_point = feature["geometry"]["coordinates"][0][:2]
        coordinates.append(first_point)

    # Cluster coordinates using KMeans algorithm
    kmeans = KMeans(n_clusters=num_clusters, random_state=0, n_init="auto").fit(
        coordinates
    )
    labels = kmeans.labels_

    # Add the cluster label property to each feature in the GeoJSON data
    for i, feature in enumerate(feature_collection["features"]):
        feature["properties"]["cluster_label"] = int(labels[i])

    # Add the colour property to each feature in the GeoJSON data
    # This should assign a different colour from the colours list to each feature within a cluster
    for i in range(num_clusters):
        colour_idx = 0
        for feature in feature_collection["features"]:
            if feature["properties"]["cluster_label"] == i:
                feature["properties"]["colour"] = colours[colour_idx % len(colours)]
                colour_idx += 1

    # Write the feature collection to the GeoJSON file
    with open(
        main_geojson, "w", encoding="utf-8"
    ) as f:
        geojson.dump(feature_collection, f, indent=4)

    # Remove duplicates from 'years' list
    years = list(set(years))
    years.sort(reverse=True)
    # print(years)

    print("Data created")

    return years, feature_collection

def save_tracks_as_elevation_profiles(feature_collection):
    """Create elevation profiles for each track in the feature collection and save them as images."""
    # Create the output directory if it doesn't exist
    output_dir = os.path.join(os.getcwd(), "elevation_profiles")
    os.makedirs(output_dir, exist_ok=True)

    # Find the maximum elevation found in all tracks
    max_elevation = 0
    for feature in feature_collection["features"]:
        coordinates = feature["geometry"]["coordinates"]
        for i in coordinates:
            if i[2] > max_elevation:
                max_elevation = i[2]

    # Round up max elevation to the nearest 100
    remainder = max_elevation % 100
    if remainder != 0:
        max_elevation += 100 - remainder

    # Create elevation profiles for each track
    for feature in feature_collection["features"]:
        track_name = feature["properties"]["name"].replace(" ", "").replace("@", "_")
        output_file = os.path.join(output_dir, track_name + ".png")

        # Check if the file already exists
        if os.path.exists(output_file):
            continue
        
        coordinates = feature["geometry"]["coordinates"]
        
        # Extract elevation and distance data
        elevations = [point[2] for point in coordinates]
        distances = [point[3] for point in coordinates]

        # Plotting elevation against total distance
        plt.figure(figsize=(8, 4))
        plt.plot(distances, elevations)
        plt.xlabel('Distance (miles)')
        plt.ylabel('Elevation (m)')
        plt.grid(True)

        # Set y-axis limits
        plt.ylim(bottom=0, top=max_elevation)

        # Save the graph as an image
        plt.savefig(output_dir + "/" + track_name)
        plt.close()  # Close the figure to prevent overlap
    
    print("Elevation profiles created")

def save_tracks_as_gpx(feature_collection):
    # Create the output directory if it doesn't exist
    output_dir = os.path.join(os.getcwd(), "track_downloads")
    os.makedirs(output_dir, exist_ok=True)

    for feature in feature_collection["features"]:
        track_name = feature["properties"]["name"].replace(" ", "").replace("@", "_")
        gpx_data = feature["geometry"]["coordinates"]

        # Create a new GPX XML document
        gpx = Element("gpx", attrib={"version": "1.1", "creator": "Your Creator Name"})
        trk = SubElement(gpx, "trk")
        trkseg = SubElement(trk, "trkseg")

        for point in gpx_data:
            trkpt = SubElement(
                trkseg, "trkpt", attrib={"lat": str(point[1]), "lon": str(point[0])}
            )
            ele = SubElement(trkpt, "ele")
            ele.text = str(point[2])  # Elevation
            # You can add more track point attributes here if needed

        # Save the GPX file
        gpx_filename = os.path.join(output_dir, f"{track_name}.gpx")
        with open(gpx_filename, "w", encoding="utf-8") as f:
            f.write(tostring(gpx, encoding="unicode"))

    print("GPX files created")

def split_features_to_files(features, output_dir):
    manifest = []
    marker_features = []
    for feature in features:
        name = feature.get("properties", {}).get("name")
        if not name:
            continue  # Skip features without a name

        # Sanitize filename:
        # Replace '@' with '_'
        safe_name = name.replace('@', '_')
        # Remove any characters that are not alphanumeric, underscore, or hyphen
        safe_name = "".join(c for c in safe_name if c.isalnum() or c in ('_', '-')).rstrip()
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
    return manifest, marker_features


def write_manifest(manifest, manifest_path):
    with open(manifest_path, "w", encoding="utf-8") as mf:
        json.dump(manifest, mf, ensure_ascii=False, indent=2)


def write_track_markers(marker_features, output_path):
    track_markers_geojson = {
        "type": "FeatureCollection",
        "features": marker_features
    }
    with open(output_path, "w", encoding="utf-8") as mf:
        json.dump(track_markers_geojson, mf, ensure_ascii=False, indent=2)

def main():
    # Paths
    main_geojson = "moorwalkers.geojson"
    output_dir = "tracks"
    manifest_path = "tracks_manifest.json"
    track_markers_path = "track_markers.geojson"

    years, feature_collection = create_data(main_geojson)

    # Create individual elevation profile images
    save_tracks_as_elevation_profiles(feature_collection)

    # Create individual gpx files from the created data for users to download
    save_tracks_as_gpx(feature_collection)

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Load GeoJSON
    with open(main_geojson, "r", encoding="utf-8") as f:
        data = json.load(f)

    features = data.get("features", [])

    manifest, marker_features = split_features_to_files(features, output_dir)
    write_manifest(manifest, manifest_path)
    write_track_markers(marker_features, track_markers_path)

    print(f"Split {len(features)} features into '{output_dir}' folder, created manifest '{manifest_path}', and created '{track_markers_path}' with {len(marker_features)} markers.")


if __name__ == "__main__":
    main()
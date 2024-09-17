'''
Sources:
Initial:
https://medium.com/@mahyar.aboutalebi/downloading-sentinel-2-imagery-in-python-with-google-colab-updated-nov-2023-f21d75a92407
Expansion for spectral bands:
https://towardsdatascience.com/satellites-can-see-invisible-lava-flows-and-active-wildires-but-how-python-371915464d1c
Sentinel 2 spectral bands:
Band 2 (Blue): 496.6 nm (10m resolution)
Band 3 (Green): 560.0 nm (10m resolution)
Band 4 (Red): 664.5 nm (10m resolution)
Band 8 (NIR — Near Infrared): 835.1 nm (10m resolution)
Band 11 (SWIR1- Shortwave Infrared): 1613.7 nm (20m resolution)
Band 12 (SWIR2- Shortwave Infrared): 2202.4 nm (20m resolution)
'''
import os
import math
import rasterio
import requests
import json
import xml.etree.ElementTree as ET
import pandas as pd
from pathlib import Path
from scipy.ndimage import zoom


def get_token():
    """
    Retrieves an access token from the Copernicus Open Access Hub for authentication.

    The function retrieves the username and password from environment variables,
    sends a POST request to the authentication server, and extracts the access token from the response.

    Parameters:
    None

    Returns:
    str: The access token obtained from the authentication server.
    """
    username = os.getenv("COPERNICUS_USERNAME")
    password = os.getenv("COPERNICUS_PASSWORD")

    auth_server_url = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
    data = {
        "client_id": "cdse-public",
        "grant_type": "password",
        "username": username,
        "password": password
    }

    response = requests.post(auth_server_url, data=data, verify=True, allow_redirects=False)
    access_token = json.loads(response.text)["access_token"]

    return access_token


def get_collections(url, satellite, level, aoi_point, start_date, end_date):
    """
    Retrieves a list of Sentinel-2 products based on specified criteria.

    Parameters:
    satellite (str): The name of the satellite (e.g., 'SENTINEL-2').
    level (str): The level of the product (e.g., 'S2MSI1C').
    aoi_point (str): The point of interest in the format "POINT(longitude latitude)".
    start_date (str): The start date for filtering products in the format "YYYY-MM-DDTHH:MM:SS.sssZ".
    end_date (str): The end date for filtering products in the format "YYYY-MM-DDTHH:MM:SS.sssZ".

    Returns:
    pandas.DataFrame: A DataFrame containing the filtered products.
    The DataFrame includes columns such as 'Id', 'Name', 'ContentDate', etc.
    """
    query = f"{url}/Products?$filter=Collection/Name eq '{satellite}' and Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'productType' and att/OData.CSC.StringAttribute/Value eq '{level}') and OData.CSC.Intersects(area=geography'SRID=4326;{aoi_point}') and ContentDate/Start gt {start_date} and ContentDate/Start lt {end_date}"
    response = requests.get(query).json()
    result = pd.DataFrame.from_dict(response["value"])

    # Filter records where 'online' column is True
    result = result[result['Online'] is True]


def get_mtd_xml(url, id, name, session, path):
    """
    Downloads the metadata XML file (MTD_MSIL1C.xml) from a Sentinel-2 product from the Copernicus Open Access Hub.

    Parameters:
    url (str): The base URL of the Copernicus Open Access Hub.
    id (str): The ID of the Sentinel-2 product.
    name (str): The name of the Sentinel-2 product.
    session (requests.Session): A session object for making HTTP requests.
    path (str): The path where the downloaded XML file will be saved.

    Returns:
    str: The path to the downloaded XML file.
    """
    url_MTD = f"{url}/Products({id})/Nodes({name})/Nodes(MTD_MSIL1C.xml)/$value"
    response = session.get(url_MTD, allow_redirects=False)
    while response.status_code in (301, 302, 303, 307):
        url_MTD_location = response.headers["Location"]
        response = session.get(url_MTD_location, allow_redirects=False)

    file = session.get(url_MTD_location, verify=False, allow_redirects=True)

    # Save the product in home directory
    outfile = Path(os.path.join(path, "content", "MTD_MSIL1C.xml"))
    outfile.write_bytes(file.content)

    return str(outfile)


def download_band(band_node, url, id, name, session, path):
    """
    Downloads a specific band of a Sentinel-2 product from the Copernicus Open Access Hub.

    Parameters:
    band_node (list): A list containing the nodes to access the band in the Sentinel-2 product.
    url (str): The base URL of the Copernicus Open Access Hub.
    id (str): The ID of the Sentinel-2 product.
    name (str): The name of the Sentinel-2 product.
    session (requests.Session): A session object for making HTTP requests.
    path (str): The path where the downloaded band will be saved.

    Returns:
    None. The function saves the downloaded band to the specified path.
    """
    url_full = f"{url}/Products({id})/Nodes({name})/Nodes({band_node[1]})/Nodes({band_node[2]})/Nodes({band_node[3]})/Nodes(R10m)/Nodes({band_node[4]})/$value"
    response = session.get(url_full, allow_redirects=True)
    while response.status_code in (301, 302, 303, 307):
        url_full = response.headers["Location"]
        response = session.get(url_full, allow_redirects=False)
    file = session.get(url_full, verify=False, allow_redirects=True)
    # Save the product
    outfile = Path(os.path.join(path, "content", band_node[4]))
    outfile.write_bytes(file.content)
    print("Saved:", band_node[4])


def calculate_new_coordinates(center_lat, center_lon, bearing, distance=3):
    """
    Calculates the new coordinates (longitude and latitude) given a center point, bearing, and distance.

    Parameters:
    center_lat (float): The latitude of the center point in degrees.
    center_lon (float): The longitude of the center point in degrees.
    bearing (float): The bearing in radians from the center point.
    distance (float, optional): The distance in kilometers from the center point. Defaults to 3 km.

    Returns:
    tuple: A tuple containing the new longitude and latitude in degrees.
    """
    # Earth radius in kilometers
    earth_radius = 6371.0

    # Convert coordinates to radians
    center_lat_rad = math.radians(center_lat)
    center_lon_rad = math.radians(center_lon)

    # Calculate new latitude
    new_lat_rad = math.asin(math.sin(center_lat_rad) * math.cos(distance / earth_radius) +
                            math.cos(center_lat_rad) * math.sin(distance / earth_radius) * math.cos(bearing))

    # Calculate new longitude
    new_lon_rad = center_lon_rad + math.atan2(math.sin(bearing) * math.sin(distance / earth_radius) *
                                              math.cos(center_lat_rad),
                                              math.cos(distance / earth_radius) - math.sin(center_lat_rad) *
                                              math.sin(new_lat_rad))

    # Convert back to degrees
    new_lat = math.degrees(new_lat_rad)
    new_lon = math.degrees(new_lon_rad)

    return new_lon, new_lat


def downscale_raster(root_path, jp2_root, band):
    """
    Resamples a raster image to a lower resolution using bilinear interpolation.

    Parameters:
    root_path (str): The root directory path where the input and output files are located.
    jp2_root (str): The root name of the input raster file (without the band suffix).
    band (str): The band name of the input raster file (e.g., 'B11', 'B12').

    Returns:
    None. The function writes the resampled raster to a new file with "_resampled" appended to the band name.
    """
    input_band_path = os.path.join(root_path, "content", f"{jp2_root}_{band}.jp2")
    output_band_path = input_band_path.replace(".jp2", "_resampled.jp2")

    scale_factor = 1/2
    with rasterio.open(input_band_path) as src:
        # Read the data
        data = src.read(1)

        # Calculate the new dimensions
        new_height = int(src.height / scale_factor)
        new_width = int(src.width / scale_factor)

        # Use scipy's zoom function for resampling
        resampled_data = zoom(data, 1/scale_factor, order=3)

        # Update metadata for the new raster
        transform = src.transform * src.transform.scale(
            (src.width / resampled_data.shape[1]),
            (src.height / resampled_data.shape[0])
        )

        new_profile = src.profile
        new_profile.update({
            'driver': 'JP2OpenJPEG',
            'height': new_height,
            'width': new_width,
            'transform': transform
        })

        # Write the resampled raster to a new file
        with rasterio.open(output_band_path, 'w', **new_profile) as dst:
            dst.write(resampled_data, 1)


def main():
    """
    Main function to download and process Sentinel-2 satellite imagery data.

    Parameters:
    None

    Returns:
    None
    """
    root_path = os.path.dirname(os.path.realpath(__file__))
    url_dataspace = "https://catalogue.dataspace.copernicus.eu/odata/v1"

    # Filtering
    satellite = "SENTINEL-2"
    level = "S2MSI1C"
    access_token = get_token()

    # Lava flow
    lava_lat = -22.411503
    lava_long = 63.892295
    lava_aoi_point = f"POINT({str(lava_lat)} {str(lava_long)})"
    lava_start_date = "2024-02-07"
    lava_end_date = "2024-02-10"
    lava_jp2_root = "T27VVL_20240208T130311"

    # Fire
    # fire_lat = -119.26
    # fire_long = 37.1914
    # fire_aoi_point = f"POINT({str(fire_lat)} {str(fire_long)})"
    # fire_start_date = "2020-09-07"
    # fire_end_date = "2020-09-10"
    # fire_jp2_root = "T11SKB_20200908T183921"

    start_date_full = lava_start_date + "T00:00:00.000Z"
    end_date_full = lava_end_date + "T00:00:00.000Z"

    result = get_collections(url_dataspace, satellite, level, lava_aoi_point, start_date_full, end_date_full)

    session = requests.Session()
    session.headers["Authorization"] = f"Bearer {access_token}"

    product_id = result.iloc[0, 1]
    product_name = result.iloc[0, 2]

    outfile = get_mtd_xml(url_dataspace, product_id, product_name, session, root_path)

    # Pass the path of the xml document
    tree = ET.parse(str(outfile))
    # get the parent tag
    root = tree.getroot()

    # Get the location of individual bands in Sentinel-2 granule
    band_path = []
    band_path.append(f"{product_name}/{root[0][0][12][0][0][0].text}.jp2".split("/"))  # Blue
    band_path.append(f"{product_name}/{root[0][0][12][0][0][1].text}.jp2".split("/"))  # Green
    band_path.append(f"{product_name}/{root[0][0][12][0][0][2].text}.jp2".split("/"))  # Red
    band_path.append(f"{product_name}/{root[0][0][12][0][0][3].text}.jp2".split("/"))  # Near-infrared
    band_path.append(f"{product_name}/{root[0][0][12][0][0][11].text}.jp2".split("/"))  # Shortwave infrared-1
    band_path.append(f"{product_name}/{root[0][0][12][0][0][12].text}.jp2".split("/"))  # Shortwave infrared-2

    # Build the url
    for band_node in band_path:
        download_band(band_node, url_dataspace, product_id, product_name, session, root_path)

    # Center coordinates
    center_lat = 63.892295
    center_lon = -22.411503

    # Calculate coordinates for the four corners
    _, north_lat = calculate_new_coordinates(center_lat, center_lon, 0)
    _, south_lat = calculate_new_coordinates(center_lat, center_lon, math.pi)
    east_lon, _ = calculate_new_coordinates(center_lat, center_lon, math.pi / 2)
    west_lon, _ = calculate_new_coordinates(center_lat, center_lon, -math.pi / 2)

    # Print the coordinates in the desired format
    print(f"({west_lon:.4f} {north_lat:.4f}, {east_lon:.4f} {north_lat:.4f}, {east_lon:.4f} {south_lat:.4f}, \
          {west_lon:.4f} {south_lat:.4f}, {west_lon:.4f} {north_lat:.4f})")

    downscale_raster(root_path, lava_jp2_root, 'B11')
    downscale_raster(root_path, lava_jp2_root, 'B12')


if __name__ == "__main__":
    main()

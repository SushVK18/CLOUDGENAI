import boto3
import os
from botocore.exceptions import ClientError

def simple_download_from_s3(bucket_name, s3_key, local_path, override=False):
    """
    Download a file from S3 to a local path.
    
    :param bucket_name: Name of the S3 bucket
    :param s3_key: S3 key (path) of the file to download
    :param local_path: Local directory where the file will be saved
    :param override: If True, overwrite the existing file; if False, abort if the file exists
    """
    # Extract the filename from the S3 key
    filename = os.path.basename(s3_key)
    full_local_path = os.path.join(local_path, filename)

    # Check if the local file already exists
    if not override and os.path.isfile(full_local_path):
        print(f"File already exists at {full_local_path}. Download aborted.")
        return None  # or raise an exception if you prefer

    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(os.path.abspath(full_local_path)), exist_ok=True)
    
    # Initialize S3 client - no credentials needed, uses IAM role
    s3_client = boto3.client('s3')
    
    # Download the file
    s3_client.download_file(bucket_name, s3_key, full_local_path)
    print(f"File downloaded successfully to: {full_local_path}")
    return full_local_path

import boto3
import os
from botocore.exceptions import ClientError

def simple_upload_to_s3(bucket_name, local_file, s3_key, override=False):
    """
    Upload a file from a local path to S3.
    
    :param bucket_name: Name of the S3 bucket
    :param local_path: Local path of the file to upload
    :param s3_key: S3 key (path) where the file will be uploaded
    :param override: If True, overwrite the existing file in S3; if False, abort if the file exists
    """
    # Check if the local file exists
    if not os.path.isfile(local_file):
        raise FileNotFoundError(f"The file {local_file} does not exist.")
    
    # Initialize S3 client - no credentials needed, uses IAM role
    s3_client = boto3.client('s3')
    
    # Check if the file already exists in S3
    if not override:
        try:
            s3_client.head_object(Bucket=bucket_name, Key=s3_key)
            print(f"File already exists at s3://{bucket_name}/{s3_key}. Upload aborted.")
            return None  # or raise an exception if you prefer
        except ClientError as e:
            # If the error is 404, the file does not exist
            if e.response['Error']['Code'] != '404':
                # If it's a different error, raise it
                raise

    # Upload the file
    s3_client.upload_file(local_file, bucket_name, s3_key)
    print(f"File uploaded successfully to s3://{bucket_name}/{s3_key}")
    return s3_key

import boto3
import os
from botocore.exceptions import ClientError

def upload_folder_to_s3(bucket_name, local_folder, s3_folder):
    """
    Upload all files from a local folder to an S3 bucket without overwriting existing files.
    
    :param bucket_name: Name of the S3 bucket
    :param local_folder: Local folder containing files to upload
    :param s3_folder: S3 folder (prefix) where the files will be uploaded
    """
    # Initialize S3 client - no credentials needed, uses IAM role
    s3_client = boto3.client('s3')

    # Iterate through all files in the local folder
    for root, _, files in os.walk(local_folder):
        for file in files:
            local_file_path = os.path.join(root, file)
            # Construct the S3 key by combining the s3_folder and the file name
            s3_key = os.path.join(s3_folder, file)

            # Check if the file already exists in S3
            try:
                s3_client.head_object(Bucket=bucket_name, Key=s3_key)
                print(f"File already exists at s3://{bucket_name}/{s3_key}. Upload aborted.")
            except ClientError as e:
                # If the error is 404, the file does not exist
                if e.response['Error']['Code'] == '404':
                    # Upload the file
                    s3_client.upload_file(local_file_path, bucket_name, s3_key)
                    print(f"File uploaded successfully to s3://{bucket_name}/{s3_key}")
                else:
                    # If it's a different error, raise it
                    raise

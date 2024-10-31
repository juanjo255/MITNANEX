#!/usr/bin/env python3
import requests
import sys

if __name__ == "__main__":
    args = sys.argv
    try:
        response = requests.post("https://mitomap.org/mitomaster/websrvc.cgi", files={"file": open(args[2]),'fileType': ('', 'snvlist'),'output': ('', 'detail')})
        #print(str(response.content, 'utf-8'))
        file_mitomap = open(args[1], "w")
        file_mitomap.write(str(response.content, 'utf-8'))
        file_mitomap.close()
    except requests.exceptions.HTTPError as err:
        print("HTTP error: " + err)
    except:
        print("Error")  
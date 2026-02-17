# Python file temporarily created to convert non-specific filetype to .txt
# since Gwyddion output was not clearly specified

# Author: Berke Santos & Giuseppe Legrotaglie
# Developed with the help of Claude.AI
# Created: 16/02/2026

import os
import sys

def add_txt_to_files(folder_path):
    if not os.path.isdir(folder_path):
        print(f"Error: '{folder_path}' is not a valid directory.")
        return
    
    count = 0
    for filename in os.listdir(folder_path):
        if filename.lower().endswith('.txt'):
            continue  # Skip existing .txt files
        
        old_path = os.path.join(folder_path, filename)
        if os.path.isfile(old_path):  # Only files, not subdirs
            new_filename = filename + '.txt'
            new_path = os.path.join(folder_path, new_filename)
            os.rename(old_path, new_path)
            print(f"Renamed: {filename} -> {new_filename}")
            count += 1
    
    print(f"Processed {count} files.")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        folder = sys.argv[1]
    else:
        folder = input("Enter folder path: ").strip()
    
    add_txt_to_files(folder)

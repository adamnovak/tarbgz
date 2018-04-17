#!/usr/bin/env python2.7

"""
tarbgz.py: random access TAR files compressed with BGZF and external indexing.

Usage examples:

    tar cf example.tar demo
    bgzip example.tar
    ./tarbgz.py example.tar.gz example.tar.gz.index --index
    ./tarbgz.py example.tar.gz example.tar.gz.index --find ""
    ./tarbgz.py example.tar.gz example.tar.gz.index --find "demo"
    ./tarbgz.py example.tar.gz example.tar.gz.index --extract "demo/test.c"

"""

import argparse
import os
import sys
import logging

from Bio import bgzf
import tarfile
import pickle
import shutil

class BgzfWrapper(object):
    """
    File-like object to wrap a Bio.bgzf BgzfReader, with ordinary
    non-virtual-offset-based forward-only seek() and tell().
    
    """
    
    def __init__(self, reader, initial_offset=0):
        """
        Wrap the given BgzfReader.
        
        If the BgzfReader is not at the start of the file, initial_offset must
        be set to the offset in the compressed data corresponding to the
        virtual offset that the BgzfReader is at (to allow absolute seek
        forward from there, as is used by tarfile).
        
        """
        
        # Save the reader
        self.reader = reader
        
        # Track offset
        self.offset = initial_offset
        
    def read(self, size=-1):
        """
        Read the given number of bytes, or whatever is available.
        Note that reading whatever is available isn't actually implemented and will fail.
        
        """
        
        bytes_read = self.reader.read(size)
        self.offset += len(bytes_read)
        
        logging.info("Read {} to {} in compressed stream".format(len(bytes_read), self.offset))
        
        return bytes_read
        
    def readline(self):
        """
        Read a single line.
        """
        
        line = self.reader.readline()
        self.offset += len(line)
        
        logging.info("Read {} to {} in compressed stream".format(len(line), self.offset))
        
        return line
        
    def tell(self):
        """
        Report the current offset from the file start.
        """
        return self.offset
        
    def seek(self, offset, whence = os.SEEK_SET):
        # How far should we go?
        distance = offset
        if whence == os.SEEK_SET:
            # We are advancing to an absolute and not a relative position
            distance -= self.offset
            
        # We can't seek from the end
        assert(whence != os.SEEK_END)
        # And we can't seek backward
        assert(distance >= 0)
        
        logging.info("Seek ahead {} from {} in compressed stream".format(distance, self.offset))
        
        while distance > 0:
            dist_read = len(self.read(distance))
            distance -= dist_read
            
            if dist_read == 0:
                # We hit EOF
                break
                
class IndexEntry(object):
    """
    Entry for a file in the index.
    
    Holds virtual and real offsets, and size, but not the file's name, which is
    implicit in the entry's storage location.

    """
    
    def __init__(self, virtual_offset, real_offset, size, next_virtual_offset):
        """
        Make a new entry with the given virtual and real offsets and file size.
        The next virtual offset is also stored, allowing the real byte range in the underlying bgzip file to be calculated.
        """
        
        self.virtual_offset = virtual_offset
        self.real_offset = real_offset
        self.size = size
        self.next_virtual_offset = next_virtual_offset
        
    def __repr__(self):
        """
        Represent this entry as a string.
        """
        
        return "IndexEntry({}, {}, {}, {})".format(self.virtual_offset, self.real_offset, self.size, self.next_virtual_offset)
        
class IndexNode:
    """
    Node in the index tree.
    
    Holds a dict of children, and also possibly an IndexEntry for the directory or file it represents.
    """
    
    def __init__(self, children=None, entry=None):
        """
        Make a new IndexNode with the given children and entry (defaulting to empty).
        """
        
        self.children = {} if children is None else children
        self.entry = entry
        
    def has_child(self, name):
        """
        Return true if we have a child node with the given name, and false otherwise.
        """
        
        return self.children.has_key(name)
        
    def make_child(self, name):
        """
        Create a new child of this node under the given name.
        Overwrites any existing child.
        """
        
        self.children[name] = IndexNode()
        
    def get_child(self, name):
        """
        Get the child with the given name, or throw an error if no such child exists.
        """
        
        return self.children[name]
        
    def get_children(self):
        """
        Return an iterator over all child names.
        """
        
        return self.children.iterkeys()
        
    def __repr__(self):
        """
        Represent this node as a string.
        """
        
        return "IndexNode({}, {})".format(self.children, self.entry)
        


class Index(object):
    """
    Represents an index of an archive. Holds IndexEntry objects in a
    hierarchical filesystem.
    """            
    
    def __init__(self, root = None):
        """
        Make a new, empty index.
        """
        
        # We store the whole index as a tree of IndexNodes
        self.root = root if root is not None else IndexNode()
        
    def insert(self, path, entry):
        """
        Insert a file or directory into the index at the given path,
        with the given IndexEntry holding its metadata.
        """
        
        # Get all the parts of the path
        parts = Index.atomize_path(path)
        
        directory = self.root
        
        for part in parts:
            if not directory.has_child(part):
                # Make sure all the parent directories exist, although they may not have IndexEntries yet.
                directory.make_child(part)
            # Advance into the thing we are looking for
            directory = directory.get_child(part)
            
        # Now the "directory" is the file or directory we are indexing.
        # Set its index entry
        directory.entry = entry
        
    def readdir(self, path):
        """
        Return an iterator over all the string names of files or directories in the given directory in the index.
        
        Use "" to look in the root.
        """
        
        # Get all the parts of the path
        parts = Index.atomize_path(path)
        
        directory = self.root
        
        for part in parts:
            if not directory.has_child(part):
                raise RuntimeError("{} not found".format(part))
            # Advance into the thing we are looking for
            directory = directory.get_child(part)
            
        return directory.get_children()
        
    def get(self, path):
        """
        Return the index entry for the given path, or None if the path has no entry or does not exist.
        """
        
        # Get all the parts of the path
        parts = Index.atomize_path(path)
        
        directory = self.root
        
        for part in parts:
            if not directory.has_child(part):
                return None
            # Advance into the thing we are looking for
            directory = directory.get_child(part)
            
        return directory.entry
        
    def __repr__(self):
        """
        Represent this index as a string.
        """
        
        return "Index({})".format(self.root)
        
    
    @staticmethod            
    def atomize_path(path):
        """
        Split a path with os.path.split into a list of path components.
        """
        
        # Split off the last path component
        (head, tail) = os.path.split(path)
        
        if head != "":
            # Recurse on the head
            parts = Index.atomize_path(head)
        else:
            # This was just a zero- or one-component path
            parts = []
        
        if tail != "":
            # This was not a 0-component path
            parts.append(tail)
            
        return parts
    
    


def parse_args(args):
    """
    Takes in the command-line arguments list (args), and returns a nice argparse
    result with fields for all the options.
    
    Borrows heavily from the argparse documentation examples:
    <http://docs.python.org/library/argparse.html>
    """
    
    # Construct the parser (which is stored in parser)
    # Module docstring lives in __doc__
    # See http://python-forum.com/pythonforum/viewtopic.php?f=3&t=36847
    # And a formatter class so our examples in the docstring look good. Isn't it
    # convenient how we already wrapped it to 80 characters?
    # See http://docs.python.org/library/argparse.html#formatter-class
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("archive_path",
        help=".tar.bgz filename to operate on")
    parser.add_argument("index_path",
        help="filename of index to read or write")
    parser.add_argument("--index", action="store_true",
        help="compute the index")
    parser.add_argument("--find", type=str, default=None,
        help="find the given file or list the given directory in the index")
    parser.add_argument("--extract", type=str, default=None,
        help="extract the given file or directory from the archive")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)

def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    logging.basicConfig(level="INFO")
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    if options.index:
        # Build an index.
        
        # This will hold IndexEntries for each path in the tar file
        index = Index()
    
        # Open the archive
        archive_bgzf = bgzf.BgzfReader(options.archive_path, "rb")
        archive_wrapper = BgzfWrapper(archive_bgzf)
        
        # Track the BGZF virtual offset of the next TAR record
        bgzf_offset = archive_bgzf.tell()
        # And the tar offset
        tar_offset = archive_wrapper.tell()
        
        # Don't actually open the tar file until we get some offsets because it immediately reads the first block
        archive_tar = tarfile.TarFile(fileobj=archive_wrapper)
        
        for info in archive_tar:
            # Dump each info record
            logging.info("Got tar info for {}".format(info.name))
            logging.info("Virtual offset: {} Decompressed offset: {}".format(bgzf_offset, tar_offset))
            logging.info("tarinfo offset: {} tarinfo data offset: {} tarinfo size: {}".format(info.offset, info.offset_data, info.size))
            logging.info("tarfile offset: {}".format(archive_tar.offset))
           
            # Before we can make a record for this file, we need to know where the next file/EOF is.
           
            # archive_tar.offset holds the uncompressed file offset of the *next* entry in the tar file.
            # Seek there so we can get its bgzf virtual offset.
            archive_wrapper.seek(archive_tar.offset)
            logging.info("Manual seek to {} complete".format(archive_tar.offset))
            
            # Make an IndexEntry for this file, but which knows the virtual offset of the next file
            # This means it can define a used byte range in the bgzf file for asynchronous retrieval
            entry = IndexEntry(bgzf_offset, tar_offset, info.size, archive_bgzf.tell())
            
            # Record data in the index
            index.insert(info.name, entry)
            
            # Update offsets for next file
            bgzf_offset = archive_bgzf.tell()
            tar_offset = archive_wrapper.tell()
            
        logging.info("Index complete")
        logging.info(index)
        
        # Save the whole index to a file.
        # TODO: Allow random access without loading the whole index by saving to a cool binary format.
        pickle.dump(index, open(options.index_path, 'w'))
            
    elif options.find is not None:
        # We want to do a list of the given directory.
        
        # Load the index
        index = pickle.load(open(options.index_path, 'r'))
        
        for basename in index.readdir(options.find):
            # For each child of what we want to list
            # Compute its full path
            full_name = os.path.join(options.find, basename)
            
            # Try to determine its size, if it has an entry
            size = 0
            entry = index.get(full_name)
            if entry is not None:
                size = entry.size
            
            # TODO: figure out if things are directories.
            print("{}\t{}".format(full_name, size))
            
    elif options.extract is not None:
        # We want to extract the given file.
        
        # Load the index
        index = pickle.load(open(options.index_path, 'r'))
        
        # Find the file to extract
        entry = index.get(options.extract)
        if entry is None:
            logging.critical("{} not found in index".format(options.extract))
            return 1
        
        # Open the archive and seek.
        # TODO: Request stuff from Glacier and then go get it and trim it,
        # updating the virtual offsets to account for the blocks not present.
        archive_bgzf = bgzf.BgzfReader(options.archive_path, "rb")
        archive_bgzf.seek(entry.virtual_offset)
        archive_wrapper = BgzfWrapper(archive_bgzf, entry.real_offset)
        
        # Start tar read at that position
        archive_tar = tarfile.TarFile(fileobj=archive_wrapper)
        
        # Iterate the first tar entry
        info = next(archive_tar)
        # And get a file object for its data
        stream = archive_tar.extractfile(info)
        
        if stream is None:
            # Not a file, so logging.info nothing
            return
        
        # Dump to standard output
        shutil.copyfileobj(stream, sys.stdout)

def entrypoint():
    """
    0-argument entry point for setuptools to call.
    """
    
    # Provide main with its arguments and handle exit codes
    sys.exit(main(sys.argv))
    
if __name__ == "__main__" :
    entrypoint()
        

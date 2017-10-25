import numpy as np


class Path(object):
    def __init__(self, np_array_view, length):
        self._path = np_array_view
        self._length = length


class Paths(object):
    def __init__(self, num_paths, max_len):
        self._paths = np.zeros((num_paths, max_len), dtype=int)
        self._path_lengths = np.zeros(num_paths, dtype=int)    

    def path_length(self, idx):
        """
        Return the length of a path.

        INPUT:
        - ``idx`` -- the index of the path

        """
        return self._path_lengths[abs(idx)-1]

    def _set_path_length(self, idx, length):
        """
        Set the length of a path.

        INPUT:
        - ``idx`` -- the index of a path
        - ``length`` -- the new length

        """
        self._path_lengths[abs(idx)-1] = length

    def get_path(self, idx):
        """
        Return a path (as a view on a numpy array)

        INPUT:
        - ``idx`` -- the index of the path. If negative, a reversed 
            path is returned.

        """
        length = self.path_length(idx)
        if idx > 0:
            return self._paths[idx-1][:length]
        else:
            return self._paths[-idx-1][length-1::-1]

    def append_to_path(self, idx, appended_path):
        """
        Append to one of the paths.

        INPUT:
        - ``idx`` -- the index of the path to append to
        - ``appended_path`` -- the path appended to the path of index
            ``idx``

        """
        length = self.path_length(idx)
        app_length = len(appended_path)
        if idx > 0:
            self._paths[idx-1][length:length+app_length] =\
                appended_path
        else:
            path = self._paths[-idx-1]
            path[app_length:app_length+length] = \
                path[:length]
            path[:app_length] = appended_path[::-1]
        self._path_lengths[abs(idx)-1] += app_length

    def replace_interval(self, path_idx, begin_pos, end_pos, inserted_path):
        """
        Replace part of a path with another path.

        INPUT:
        - ``path_idx`` -- the index of the path
        - ``begin_pos`` -- the starting position of the interval (inclusive)
        - ``end_pos`` -- the ending position of the interval (exclusive)
        - ``inserted_path`` -- list or numpy array to be inserted in place of the interval ``[begin_pos:end_pos]`` of our path

        EXAMPLES:
        >>> from macaw.train_tracks.dehn_thurston.paths import Paths
        >>> paths = Paths(5, 10)
        >>> paths.replace_interval(3, 0, 0, [1, 2, 3, 4, 5, 6])
        >>> all(paths.get_path(3) == [1, 2, 3, 4, 5, 6])
        True
        >>> paths.replace_interval(3, 2, 4, [100, 101, 102])
        >>> all(paths.get_path(3) == [1, 2, 100, 101, 102, 5, 6])
        True
        >>> paths.replace_interval(-3, 0, 6, [204])
        >>> all(paths.get_path(3) == [1, 204])
        True
        >>> paths.replace_interval(-3, 1, 2, [])
        >>> all(paths.get_path(3) == [204])
        True

        """
        length = self.path_length(path_idx)
        ins_length = len(inserted_path)
        tail_length = length-end_pos
        new_length = begin_pos+ins_length+tail_length
        path = self._paths[abs(path_idx)-1]
        if path_idx > 0:
            path[begin_pos+ins_length:new_length] =\
                path[end_pos:length]
            path[begin_pos:begin_pos+ins_length] = \
                inserted_path
        else:
            path[new_length-begin_pos:new_length] =\
                path[length-begin_pos:length]
            path[length-end_pos:length-end_pos+ins_length] =\
                inserted_path[::-1]
        self._path_lengths[abs(path_idx)-1] += ins_length-(end_pos-begin_pos)



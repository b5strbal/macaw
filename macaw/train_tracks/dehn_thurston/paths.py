import numpy as np


class Path(object):
    def __init__(self, np_array_view, length_view, is_reversed):
        self._path_view = np_array_view
        self._length_view = length_view
        self._is_reversed = is_reversed

    def length(self):
        """
        Return the length of the path.
        """
        return self._length_view[0]

    def view(self):
        """
        Return a view on the path.
        """
        temp = self._path_view[:self.length()]
        if self._is_reversed:
            return temp[::-1]
        else:
            return temp

    def reverse(self):
        """
        Reverse self.
        """
        self._is_reversed = not self._is_reversed

    def replace_interval(self, begin_pos, end_pos, inserted_path):
        """
        Replace part of the path with a different sequence. 

        INPUT:
        - ``begin_pos`` -- the starting position of the interval (inclusive)
        - ``end_pos`` -- the ending position of the interval (exclusive)
        - ``inserted_path`` -- list or numpy array to be inserted in place of the interval ``[begin_pos:end_pos]`` of our path

        EXAMPLES:
        >>> from macaw.train_tracks.dehn_thurston.paths import Path
        >>> import numpy as np
        >>> arr = np.zeros(10, dtype=int)
        >>> length_arr = np.zeros(1, dtype=int)
        >>> path = Path(arr, length_arr, False)
        >>> path.replace_interval(0, 0, [1, 2, 3, 4, 5, 6])
        >>> all(path.view() == [1, 2, 3, 4, 5, 6])
        True
        >>> path.replace_interval(2, 4, [100, 101, 102])
        >>> all(path.view() == [1, 2, 100, 101, 102, 5, 6])
        True
        >>> path.reverse()
        >>> all(path.view() == [6, 5, 102, 101, 100, 2, 1])
        True
        >>> path.replace_interval(0, 6, [204])
        >>> all(path.view() == [204, 1])
        True
        >>> path.replace_interval(1, 2, [])
        >>> all(path.view() == [204])
        True

        """
        length = self.length()
        ins_length = len(inserted_path)
        tail_length = length-end_pos
        new_length = begin_pos+ins_length+tail_length
        path = self._path_view
        if not self._is_reversed:
            path[begin_pos+ins_length:new_length] =\
                path[end_pos:length]
            path[begin_pos:begin_pos+ins_length] = \
                inserted_path
        else:
            path[new_length-begin_pos:new_length] =\
                path[length-begin_pos:length]
            path[length-end_pos:length-end_pos+ins_length] =\
                inserted_path[::-1]
        self._length_view[0] += ins_length-(end_pos-begin_pos)

    def append_to_path(self, appended_path):
        """
        Append to a sequence to the path.

        INPUT:
        - ``appended_path`` -- the sequence to append

        """
        self.replace_interval(self.length(), self.length(), appended_path)

    def delete_from_path(self, begin_pos, end_pos):
        """
        Deletes part of the path.

        INPUT:
        - ``begin_pos`` -- the starting position of the interval (inclusive)
        - ``end_pos`` -- the ending position of the interval (exclusive)

        """
        self.replace_interval(begin_pos, end_pos, [])

class Paths(object):
    def __init__(self, num_paths, max_len):
        self._paths = np.zeros((num_paths, max_len), dtype=int)
        self._path_lengths = np.zeros(num_paths, dtype=int)    

    # def get_length(self, idx):
    #     """
    #     Return the length of the path.

    #     INPUT:
    #     - ``idx`` -- the index of the path

    #     """
    #     return self._path_lengths[idx-1]

    # def _set_path_length(self, idx, length):
    #     """
    #     Set the length of a path.

    #     INPUT:
    #     - ``idx`` -- the index of a path
    #     - ``length`` -- the new length

    #     """
    #     self._path_view_lengths[abs(idx)-1] = length

    def get_path(self, idx):
        """
        Return a path (as a view on a numpy array)

        INPUT:
        - ``idx`` -- the index of the path. If negative, a reversed 
            path is returned.

        """
        if idx > 0:
            is_reversed = False
        elif idx < 0:
            is_reversed = True
        else:
            assert False
        return Path(self._paths[abs(idx)-1], 
            self._path_lengths[abs(idx)-1:abs(idx)],
            is_reversed)


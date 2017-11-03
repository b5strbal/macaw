import numpy as np


class NumpyList(object):
    """
    A list-like container implemented using numpy.
    """
    def __init__(self, np_array, length_view, is_reversed):
        """
        INPUT:
        - ``np_array`` -- a 1-dimensional numpy array, containing the elements of the list
        - ``length_view`` -- a 1-dimensional numpy array of length 1 containing the length of the list. Only the first ``length_view[0]`` elements of ``np_array`` are considered to be part of the list.
        - ``is_reversed`` -- True or False. If True, the list is considered to start at the element of ``np_array`` of index ``length_view[0]-1`` and end at the element of index 0. If False, then ``np_array`` is read in the usual order.
        """
        self._array = np_array
        self._length_view = length_view
        self._is_reversed = is_reversed

    def __iter__(self):
        return self

    def __getitem__(self, key):
        return self._array[key]

    def __setitem__(self, key, item):
        self._array[key] = item

    def __repr__(self):
        return repr(self.view())

    def length(self):
        """
        Return the length of the list.

        EXAMPLES:
        >>> from macaw.train_tracks.dehn_thurston.numpy_list import NumpyList
        >>> import numpy as np
        >>> mylist = NumpyList(np.array([1,2,0,0]), np.array([3]), False)
        >>> mylist.length()
        3

        """
        return self._length_view[0]

    def view(self):
        """
        Return a view on the list.

        OUTPUT:
        a 1D numpy array containing the elements of our list. This array should always be read from left to right.

        EXAMPLES:
        >>> from macaw.train_tracks.dehn_thurston.numpy_list import NumpyList
        >>> import numpy as np
        >>> mylist = NumpyList(np.array([1,2,0,0]), np.array([3]), False)
        >>> all(mylist.view() == [1, 2, 0])
        True
        >>> mylist = NumpyList(np.array([1,2,0,0]), np.array([3]), True)
        >>> all(mylist.view() == [0, 2, 1])
        True

        """
        temp = self._array[:self.length()]
        if self._is_reversed:
            return temp[::-1]
        else:
            return temp

    def reverse(self):
        """
        Reverse self.

        EXAMPLES:
        >>> from macaw.train_tracks.dehn_thurston.numpy_list import NumpyList
        >>> import numpy as np
        >>> mylist = NumpyList(np.array([1,2,3,0]), np.array([3]), False)
        >>> all(mylist.view() == [1, 2, 3])
        True
        >>> mylist.reverse()
        >>> all(mylist.view() == [3, 2, 1])
        True

        """
        self._is_reversed = not self._is_reversed

    def replace_interval(self, begin_pos, end_pos, inserted_list):
        """
        Replace part of ``self`` with a different list.

        INPUT:
        - ``begin_pos`` -- the starting position of the interval (inclusive)
        - ``end_pos`` -- the ending position of the interval (exclusive)
        - ``inserted_list`` -- list or numpy array to be inserted in place of the interval ``[begin_pos:end_pos]`` of our list

        EXAMPLES:
        >>> from macaw.train_tracks.dehn_thurston.numpy_list import NumpyList
        >>> import numpy as np
        >>> arr = np.zeros(10, dtype=int)
        >>> length_arr = np.zeros(1, dtype=int)
        >>> mylist = NumpyList(arr, length_arr, False)
        >>> mylist.replace_interval(0, 0, [1, 2, 3, 4, 5, 6])
        >>> all(mylist.view() == [1, 2, 3, 4, 5, 6])
        True
        >>> mylist.replace_interval(2, 4, [100, 101, 102])
        >>> all(mylist.view() == [1, 2, 100, 101, 102, 5, 6])
        True
        >>> mylist.reverse()
        >>> all(mylist.view() == [6, 5, 102, 101, 100, 2, 1])
        True
        >>> mylist.replace_interval(0, 6, [204])
        >>> all(mylist.view() == [204, 1])
        True
        >>> mylist.replace_interval(1, 2, [])
        >>> all(mylist.view() == [204])
        True

        """
        length = self.length()
        ins_length = len(inserted_list)
        tail_length = length-end_pos
        new_length = begin_pos+ins_length+tail_length
        array = self._array
        if not self._is_reversed:
            array[begin_pos+ins_length:new_length] =\
                array[end_pos:length]
            array[begin_pos:begin_pos+ins_length] = \
                inserted_list
        else:
            array[new_length-begin_pos:new_length] =\
                array[length-begin_pos:length]
            array[length-end_pos:length-end_pos+ins_length] =\
                inserted_list[::-1]
        self._length_view[0] += ins_length-(end_pos-begin_pos)

    def append(self, appended_list):
        """
        Append a list to ``self``.

        INPUT:
        - ``appended_list`` -- the list to append

        """
        self.replace_interval(self.length(), self.length(), appended_list)

    def delete_interval(self, begin_pos, end_pos):
        """
        Deletes a part of the list in the specified interval

        INPUT:
        - ``begin_pos`` -- the starting position of the interval (inclusive)
        - ``end_pos`` -- the ending position of the interval (exclusive)

        """
        self.replace_interval(begin_pos, end_pos, [])


class ManyLists(object):
    """
    A collection of NumpyLists.

    The NumpyLists are index by consecutive integers from 1
    """
    def __init__(self, num_lists, max_len):
        """
        Initializes a collection of empty NumpyLists.

        INPUT:
        - ``num_lists`` -- the number of lists contained. The lists are indexed from 1 to ``num_lists``. (Not from 0!)
        - ``max_len`` -- the maximum length of the lists. Since internally the data is stored in a 2D numpy array, the maximum length for all lists is the same.

        """
        self._lists = np.zeros((num_lists, max_len), dtype=int)
        self._list_lengths = np.zeros(num_lists, dtype=int)

    def get_list(self, idx):
        """
        Return the NumpyList at the specified index.

        INPUT:
        - ``idx`` -- the index of the list. If negative, a reversed
            list is returned.

        """
        if idx > 0:
            is_reversed = False
        elif idx < 0:
            is_reversed = True
        else:
            assert False
        return NumpyList(self._lists[abs(idx)-1],
            self._list_lengths[abs(idx)-1:abs(idx)],
            is_reversed)

"""
Unit tests for plot utility functions.
"""

import unittest
import sys
import os
from unittest.mock import MagicMock

# Add project root to sys.path to allow importing 'utils'
sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), '../../'))
)

from utils.plot_utils import abbreviate_labels, apply_axis_style


class TestPlotUtils(unittest.TestCase):
    """Test cases for plotting utility functions."""

    def test_abbreviate_labels_basic(self):
        """Test standard abbreviation logic."""
        label1 = MagicMock()
        label1.get_text.return_value = "LongLabelName"
        label2 = MagicMock()
        label2.get_text.return_value = "Short"

        result = abbreviate_labels([label1, label2], 4)
        expected = ["Long", "Shor"]
        self.assertEqual(result, expected)

    def test_abbreviate_empty_list(self):
        """Test that an empty list returns an empty list (no crash)."""
        result = abbreviate_labels([], 5)
        self.assertEqual(result, [])

    def test_abbreviate_limit_zero(self):
        """Test that limit 0 returns empty strings."""
        label = MagicMock()
        label.get_text.return_value = "Anything"

        result = abbreviate_labels([label], 0)
        self.assertEqual(result, [""])

    def test_abbreviate_unicode(self):
        """Test that unicode/emojis are handled correctly."""
        label = MagicMock()
        label.get_text.return_value = "Gene_α_β"

        # Python handles unicode chars as length 1
        result = abbreviate_labels([label], 6)
        self.assertEqual(result, ["Gene_α"])

    def test_apply_axis_style(self):
        """Test that style settings are applied to labels."""
        label1 = MagicMock()
        labels = [label1]

        apply_axis_style(labels, 12, "Arial")

        label1.set_fontsize.assert_called_with(12)
        label1.set_fontfamily.assert_called_with("Arial")

    def test_apply_axis_style_empty(self):
        """Test that passing an empty list doesn't crash."""
        try:
            apply_axis_style([], 10)
        except Exception as e:
            self.fail(f"apply_axis_style crashed on empty list: {e}")


if __name__ == "__main__":
    unittest.main()
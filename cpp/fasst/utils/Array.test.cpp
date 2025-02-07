#include "Array.h"

#include <gtest/gtest.h>

using namespace fasst::utils;

// Test for fixed-size array.
TEST(FixedMatrixTest, FillAndAccess) {
    // ...existing code...
    Array<int, 3, 4> m;
    m.fill(7);
    EXPECT_EQ(m.rows(), 3);
    EXPECT_EQ(m.cols(), 4);
    EXPECT_EQ(m.size(), 12);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            EXPECT_EQ(m(i, j), 7);
        }
    }
}

TEST(FixedMatrixTest, ZeroMatrix) {
    // ...existing code...
    auto m = Array<int, 3, 2>::Zero();
    EXPECT_EQ(m.rows(), 3);
    EXPECT_EQ(m.cols(), 2);
    EXPECT_EQ(m.size(), 6);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 2; ++j) {
            EXPECT_EQ(m(i, j), 0);
        }
    }
}

TEST(FixedMatrixTest, Access) {
    // ...existing code...
    Array<int, 2, 3> m;
    m(0, 0) = 1;
    m(0, 1) = 2;
    m(0, 2) = 3;
    m(1, 0) = 4;
    m(1, 1) = 5;
    m(1, 2) = 6;
    EXPECT_EQ(m(0, 0), 1);
    EXPECT_EQ(m(0, 1), 2);
    EXPECT_EQ(m(0, 2), 3);
    EXPECT_EQ(m(1, 0), 4);
    EXPECT_EQ(m(1, 1), 5);
    EXPECT_EQ(m(1, 2), 6);
}

TEST(FixedArrayTest, AggregateInitialization) {
    // Using aggregate initialization for fixed-size Array:
    // Initialize the public data member directly.
    fasst::utils::Array<int, 2, 3> arr = { {1, 2, 3, 4, 5, 6} };
    EXPECT_EQ(arr.rows(), 2);
    EXPECT_EQ(arr.cols(), 3);
    EXPECT_EQ(arr.size(), 6);
    EXPECT_EQ(arr(0, 0), 1);
    EXPECT_EQ(arr(0, 1), 2);
    EXPECT_EQ(arr(0, 2), 3);
    EXPECT_EQ(arr(1, 0), 4);
    EXPECT_EQ(arr(1, 1), 5);
    EXPECT_EQ(arr(1, 2), 6);
}

TEST(FixedArrayTest, StringArray) {
    // Using aggregate initialization for fixed-size Array:
    // Initialize the public data member directly.
    fasst::utils::Array<std::string, 2, 3> arr = { {"one", "two", "three", "four", "five", "six"} };
    EXPECT_EQ(arr.rows(), 2);
    EXPECT_EQ(arr.cols(), 3);
    EXPECT_EQ(arr.size(), 6);
    EXPECT_EQ(arr(0, 0), "one");
    EXPECT_EQ(arr(0, 1), "two");
    EXPECT_EQ(arr(0, 2), "three");
    EXPECT_EQ(arr(1, 0), "four");
    EXPECT_EQ(arr(1, 1), "five");
    EXPECT_EQ(arr(1, 2), "six");
}

TEST(DynamicArrayTest, InitializerListInitialization) {
    // Using the initializer-list constructor for dynamic Array.
    fasst::utils::Array<int, -1, -1> dynArr(2, 3, {10, 20, 30, 40, 50, 60});
    EXPECT_EQ(dynArr.rows(), 2);
    EXPECT_EQ(dynArr.cols(), 3);
    EXPECT_EQ(dynArr.size(), 6);
    EXPECT_EQ(dynArr(0, 0), 10);
    EXPECT_EQ(dynArr(0, 1), 20);
    EXPECT_EQ(dynArr(0, 2), 30);
    EXPECT_EQ(dynArr(1, 0), 40);
    EXPECT_EQ(dynArr(1, 1), 50);
    EXPECT_EQ(dynArr(1, 2), 60);
}

// Test for dynamic array (both dimensions dynamic: use -1).
TEST(DynamicMatrixTest, ResizeAndFill) {
    // ...existing code...
    Array<int, -1, -1> m;
    m.resize(4, 4);
    m.fill(9);
    EXPECT_EQ(m.rows(), 4);
    EXPECT_EQ(m.cols(), 4);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            EXPECT_EQ(m(i, j), 9);
        }
    }
}

TEST(DynamicMatrixTest, ZeroMatrix) {
    // ...existing code...
    auto m = Array<int, -1, -1>::Zero(3, 5);
    EXPECT_EQ(m.rows(), 3);
    EXPECT_EQ(m.cols(), 5);
    EXPECT_EQ(m.size(), 15);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 5; ++j) {
            EXPECT_EQ(m(i, j), 0);
        }
    }
}

TEST(VectorTest, FixedVector) {
    // ...existing code...
    ArrayX<int, 5> v;
    v.fill(3);
    EXPECT_EQ(v.rows(), 5);
    EXPECT_EQ(v.cols(), 1);
    EXPECT_EQ(v.size(), 5);
    for (int i = 0; i < v.size(); ++i) {
        EXPECT_EQ(v(i), 3);
    }
}

TEST(VectorTest, DynamicVector) {
    // ...existing code...
    ArrayX<int, -1> v(4);
    v.fill(8);
    EXPECT_EQ(v.rows(), 4);
    EXPECT_EQ(v.cols(), 1);
    EXPECT_EQ(v.size(), 4);
    for (int i = 0; i < v.size(); ++i) {
        EXPECT_EQ(v(i), 8);
    }
    // Test resizing for dynamic vector
    v.resize(6);
    v.fill(5);
    EXPECT_EQ(v.rows(), 6);
    for (int i = 0; i < v.size(); ++i) {
        EXPECT_EQ(v(i), 5);
    }
}

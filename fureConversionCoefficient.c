#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#define U64 uint64_t
#define U32 uint32_t
#define U16 uint16_t
#define U8 uint8_t
#define loop(i, a, b) for (U32 i = a; i < b; i++)
#define glo_n 30
// const U32 glo_length = 1 << glo_n;
#define glo_length 1073741824ULL // 20 - 1048576 // 30 1073741824  //32 4294967296
void supportZhegalkin(U8 *lst_coefficients, U8 idx1)
{
    if (idx1 == glo_n)
    {
        return;
    }
    U32 half_len = 1 << idx1;
    U32 full_len = half_len << 1;
    U32 num_turns = glo_length / full_len;
    loop(idx2, 0, num_turns)
    {
        loop(idx3, 0, half_len)
            lst_coefficients[idx2 * full_len + half_len + idx3] ^= lst_coefficients[idx2 * full_len + idx3];
    }
    // loop(i, 0, glo_length)
    //     printf("%d ", lst_coefficients[i]);
    // printf("\n");
    supportZhegalkin(lst_coefficients, idx1 + 1);
}

void fastFindCoefficientsZhegalkin(U8 *lst_coefficients)
{
    U32 half_len;
    U32 full_len;
    U32 num_turns;
    loop(idx1, 0, glo_n)
    {
        half_len = 1 << idx1;
        full_len = half_len << 1;
        num_turns = glo_length / full_len;
        loop(idx2, 0, num_turns)
        {
            loop(idx3, 0, half_len)
                lst_coefficients[idx2 * full_len + half_len + idx3] ^= lst_coefficients[idx2 * full_len + idx3];
        }
    }
}

void fastFindCoefficientsReal(int8_t *lst_coefficients)
{
    U32 half_len;
    U32 full_len;
    U32 num_turns;
    loop(idx1, 0, glo_n)
    {
        half_len = 1 << idx1;
        full_len = half_len << 1;
        num_turns = glo_length / full_len;
        loop(idx2, 0, num_turns)
        {
            loop(idx3, 0, half_len)
                lst_coefficients[idx2 * full_len + half_len + idx3] -= lst_coefficients[idx2 * full_len + idx3];
        }
    }
}

void fastFindCoefficientsFourier(int8_t *lst_coefficients)
{
    U32 half_len;
    U32 full_len;
    U32 num_turns;
    loop(idx1, 0, glo_n)
    {
        half_len = 1 << idx1;
        full_len = half_len << 1;
        num_turns = glo_length / full_len;
        loop(idx2, 0, num_turns)
        {
            loop(idx3, 0, half_len)
            {
                lst_coefficients[idx2 * full_len + idx3] += lst_coefficients[idx2 * full_len + idx3 + half_len];
                lst_coefficients[idx2 * full_len + half_len + idx3] = lst_coefficients[idx2 * full_len + idx3] - 2 * lst_coefficients[idx2 * full_len + half_len + idx3];
            }
        }
    }
}
U8 vesFunction(U32 f)
{
    U8 count = 0;
    while (f)
    {
        count += f & 0b1;
        f >>= 1;
    }
    return count;
}
U8 *calcF1()
{
    U8 *table = (U8 *)malloc(sizeof(U8) * glo_length);
    U8 vec = 0;
    U8 hafl_length = glo_n >> 1;
    loop(i, 0, glo_length)
    {
        vec = vesFunction(i);
        table[i] = (vec >= hafl_length);
    }
    return table;
}

U8 *calcF2()
{
    U8 *table = (U8 *)malloc(sizeof(U8) * glo_length);
    table[0] = 1;
    table[1] = 0;
    U32 current_leng = 2;
    while (current_leng < glo_length)
    {
        loop(i, 0, current_leng)
            table[i + current_leng] = table[i] ^ 1;
        current_leng <<= 1;
    }
    return table;
}

U8 *calcF1Supper()
{
    // U8 *table = (U8 *)malloc(sizeof(U8) * glo_length);
    U8 *table = (U8 *)calloc(glo_length, sizeof(U8)); // Cấp phát bộ nhớ và khởi tạo giá trị 0 cho mảng

    U8 half_n = glo_n >> 1;
    U8 k = half_n - 1;
    U8 buff = (1 << k) - 1;
    U8 start = 1 << k;
    U8 end = 1 << (k + 1);
    while (k < glo_n)
    {
        loop(i, start + buff, end)
            table[i] = 1;
        start = end;
        end <<= 1;
        ++k;
    }
    return table;
}
U16 c = 0;
U32 f = 0;
// U8 *global_table = (U8 *)calloc(glo_length, sizeof(U8));
// Gán giá trị 1 cho tất cả các phần tử trong mảng bằng hàm memset
U8 *global_table;
void create_dynamic_array()
{
    global_table = (U8 *)malloc(glo_length * sizeof(U8)); // Cấp phát bộ nhớ cho mảng
    memset(global_table, 1, glo_length * sizeof(U8));
}
void recursive(U32 t1, U16 num_loop)
{
    if (num_loop)
    {
        U32 t2 = 1 << num_loop;
        while (t2 < t1)
        {
            f ^= t2;
            recursive(t2, num_loop - 1);
            f ^= t2;
            t2 <<= 1;
        }
    }
    else
    {
        U32 t2 = 1;
        while (t2 < t1)
        {
            f ^= t2;
            // printf("%u,", f);
            global_table[f] = 0;
            f ^= t2;
            t2 <<= 1;
            ++c;
        }
    }
}
void calcF1SpAuto()
{
    global_table[0] = 0;
    loop(t, 0, (glo_n >> 1) - 1)
    {
        c = 0;
        recursive(glo_length, t);
        // printf("t = %u ; c = %u\n", t, c);
    }
}
void calcF1Sp()
{

    // U32 limit = glo_length >> 1;
    U32 t1 = 2, t2;

    while (t1 < glo_length)
    {
        t2 = 1;
        while (t2 < t1)
        {
            printf("%u,", t2 + t1);
            ++c;
            t2 <<= 1;
        }
        t1 <<= 1;
    }
    printf("\n%u\n", c);
}
void checkF1()
{
    U8 *table1 = calcF1();
    calcF1SpAuto();
    // U8 *table2 = calcF1Supper();
    loop(i, 0, glo_length) if (table1[i] != global_table[i])
    {
        printf("failed %u %u %u\n", i, table1[i], global_table[i]);
    }
    free(table1);
    // free(table2);
}
U8 *copyArray(U8 *arr, U32 size)
{
    U8 *arrCp = (U8 *)malloc(sizeof(U8) * size);
    loop(i, 0, size)
        arrCp[i] = arr[i];
    ////
    return arrCp;
}
void printArray(int8_t *arr)
{
    loop(i, 0, glo_length)
        printf("%d,", arr[i]);
    printf("\n");
}
void printArray(U8 *arr)
{
    loop(i, 0, glo_length)
        printf("%d,", arr[i]);
    printf("\n");
}
int main()
{
    // int8_t *table = (U8 *)malloc(sizeof(int8_t) * glo_length);
    // int8_t arrTable[] = {0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1};
    //   Code mà bạn muốn đo thời gian ở đây
    clock_t start, end;
    double cpu_time_used;

    // giai đoạn 1 :
    start = clock();
    create_dynamic_array();
    calcF1SpAuto();
    // U8 *table = calcF1Supper();
    //  printArray(table);
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Thời gian thực thi giai đoạn 1: %f giây\n", cpu_time_used);
    printf("\n");

    // giai đoạn 2
    start = clock();
    fastFindCoefficientsZhegalkin(global_table);
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Thời gian thực thi giai đoạn 2: %f giây\n", cpu_time_used);
    printf("\n");
    // printArray(arrTable);
    free(global_table);
    return 0;
}
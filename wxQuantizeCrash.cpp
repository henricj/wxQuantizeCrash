
// http://trac.wxwidgets.org/ticket/17764

#include "stdafx.h"
#include "wxQuantizeCrash.h"
#include "Checksum.h"

static void Misbehave()
{
	const int ROWS = 20; // 720;
	const int COLS = 32; // 642;
	const int INPUT_SIZE = 3 * ROWS * COLS;
	const int COLORS = 16;

	JSAMPLE* input = new JSAMPLE[INPUT_SIZE];
	JSAMPLE output[ROWS * COLS];
	JSAMPROW input_row[ROWS];
	JSAMPROW output_row[ROWS];

	std::mt19937_64 engine{ ROWS * COLS };
	std::uniform_int<JSAMPLE> ui{ 0, 0xff };

	std::generate_n(input, INPUT_SIZE, [&] { return ui(engine); });

	BCryptSha256 sha256;

	unsigned char input_hash[32];

	sha256.ComputeHash(input, INPUT_SIZE, input_hash, sizeof(input_hash));

	for (auto i = 0; i < ROWS; ++i) {
		input_row[i] = &input[3 * COLS * i];
		output_row[i] = &output[COLS * i];
	}

	unsigned char palette[3 * COLORS];

	DoQuantize(COLS, ROWS, input_row, output_row, palette, COLORS);

	unsigned char input_hash2[32];

	sha256.ComputeHash(input, INPUT_SIZE, input_hash2, sizeof(input_hash2));

	if (sizeof(input_hash) != sizeof(input_hash2)
		|| 0 != memcmp(input_hash, input_hash2, sizeof(input_hash)))
	{
		std::cout << "Input changed" << std::endl;
	}

	free(input);

	const unsigned char ExpectedInputHash[32] = {
		0xeb, 0x9d, 0x67, 0x40, 0xb6, 0x69, 0xd4, 0x85, 0xfd, 0x56, 0x24, 0x6f, 0xb3, 0xa4, 0x23, 0x4a,
		0x4c, 0xfb, 0x7f, 0xbd, 0x5f, 0x53, 0x9d, 0x5b, 0x6d, 0xd6, 0x0c, 0x85, 0xe5, 0xe1, 0xac, 0x1f
	};

	const unsigned char ExpectedPaletteHash[32] = {
		0x7f, 0x25, 0xc8, 0x09, 0x1c, 0x65, 0x63, 0xce, 0x64, 0x70, 0x05, 0x62, 0xe4, 0x00, 0x31, 0xa4,
		0x3e, 0x05, 0xb3, 0xc1, 0x3d, 0x47, 0xb7, 0x9c, 0x12, 0x1e, 0xa1, 0xbb, 0x46, 0x12, 0xfb, 0x75
	};

	const unsigned char ExpectedOutputHash[32] = {
		0xf2, 0x2f, 0xb8, 0x4f, 0xde, 0xcf, 0x60, 0x17, 0x03, 0x05, 0x07, 0x1d, 0xd9, 0x6d, 0x3f, 0x0c,
		0x19, 0x95, 0x2e, 0xde, 0x1f, 0xe2, 0xa3, 0xc5, 0xd6, 0xe0, 0xbe, 0x4b, 0x88, 0xf0, 0xf2, 0xb9
	};

	if (sizeof(input_hash) != sizeof(ExpectedInputHash)
		|| 0 != memcmp(input_hash, ExpectedInputHash, sizeof(input_hash)))
	{
		std::cout << "Input mismatch" << std::endl;
	}

	unsigned char palette_hash[sizeof(ExpectedPaletteHash)];
	unsigned char output_hash[sizeof(ExpectedOutputHash)];

	sha256.ComputeHash(palette, sizeof(palette), palette_hash, sizeof(palette_hash));

	if (sizeof(ExpectedPaletteHash) != sizeof(palette_hash)
		|| 0 != memcmp(palette_hash, ExpectedPaletteHash, sizeof(ExpectedPaletteHash)))
	{
		std::cout << "Palette mismatch" << std::endl;
	}

	sha256.ComputeHash(output, sizeof(output), output_hash, sizeof(output_hash));

	if (sizeof(ExpectedOutputHash) != sizeof(output_hash)
		|| 0 != memcmp(output_hash, ExpectedOutputHash, sizeof(ExpectedOutputHash)))
	{
		std::cout << "Output mismatch" << std::endl;
	}

	{
		std::ios::fmtflags oldflags(std::cout.flags());

		std::cout << std::hex << std::setfill('0');

		for (auto y = 0; y < ROWS; ++y)
		{
			for (auto x = 0; x < COLS; ++x)
			{
				std::cout << std::setw(2) << static_cast<int>(output[y * COLS + x]) << " ";
			}

			std::cout << std::endl;
		}

		std::cout.flags(oldflags);
	}
}

int main(int argc, char *argv[]) {
	std::cout << "_MSC_FULL_VER: " << _MSC_FULL_VER
#ifdef _M_X64
		<< " x64"
#endif
#ifdef _DEBUG
		<< " Debug"
#else
		<< " Release"
#endif
		<< std::endl;

	Misbehave();

	return 0;
}

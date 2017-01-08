#include "stdafx.h"
#include "Checksum.h"

#pragma comment(lib, "bcrypt.lib")

#ifndef NT_SUCCESS
#define NT_SUCCESS(status) ((NTSTATUS)(status) >= 0)
#endif

HashHandle::HashHandle(BCRYPT_ALG_HANDLE algorithm)
{
	auto status = BCryptCreateHash(algorithm, &handle_, NULL, 0, NULL, 0, 0);
	if (!NT_SUCCESS(status))
		throw std::runtime_error("Unable to create hash");
}

HashHandle::~HashHandle()
{
	if (nullptr != handle_)
		BCryptDestroyHash(handle_);
}

BCryptAlgorithmProvider::BCryptAlgorithmProvider(LPCWSTR algorithm, LPCWSTR implementation, ULONG flags)
{
	auto ntstatus = BCryptOpenAlgorithmProvider(&handle_,
		algorithm,
		implementation,
		flags);
	if (!NT_SUCCESS(ntstatus))
		throw std::runtime_error("Unable to open algorithm provider");
}

BCryptAlgorithmProvider::~BCryptAlgorithmProvider()
{
	if (nullptr != handle_)
		BCryptCloseAlgorithmProvider(handle_, 0);
}

BCryptHash::BCryptHash(LPCWSTR algorithm, LPCWSTR implementation, ULONG flags) 
	: BCryptAlgorithmProvider(algorithm, implementation, flags)
{
	DWORD HashLength;
	DWORD ResultLength;

	auto ntstatus = BCryptGetProperty(handle(), BCRYPT_HASH_LENGTH,
		(PBYTE)&HashLength, sizeof(HashLength), &ResultLength, 0);
	if (!NT_SUCCESS(ntstatus))
		throw std::runtime_error("Unable to get hash length");

	hash_length_ = HashLength;
}

void BCryptHash::ComputeHash(const void * buffer, size_t length, void * hash, size_t hash_length)
{
	HashHandle handle{ this->handle() };

	auto status = BCryptHashData(handle.handle(), static_cast<PUCHAR>(const_cast<void*>(buffer)), static_cast<ULONG>(length), 0);
	if (!NT_SUCCESS(status))
		throw std::runtime_error("Unable to hash the data");

	status = BCryptFinishHash(handle.handle(), static_cast<PUCHAR>(hash), static_cast<ULONG>(hash_length), 0);
	if (!NT_SUCCESS(status))
		throw std::runtime_error("Unable to get the hash");
}


#pragma once

class HashHandle
{
	BCRYPT_HASH_HANDLE  handle_;
public:
	HashHandle(BCRYPT_ALG_HANDLE algorithm);
	~HashHandle();

	BCRYPT_HASH_HANDLE handle() const { return handle_; }

	HashHandle() = delete;
	HashHandle(const HashHandle&) = delete;
	HashHandle& operator=(const HashHandle&) = delete;
};

class BCryptAlgorithmProvider
{
	BCRYPT_ALG_HANDLE   handle_;
public:
	BCRYPT_ALG_HANDLE handle() const { return handle_; }

	BCryptAlgorithmProvider(LPCWSTR algorithm, LPCWSTR implementation, ULONG flags);
	~BCryptAlgorithmProvider();

	BCryptAlgorithmProvider() = delete;
	BCryptAlgorithmProvider(const BCryptAlgorithmProvider&) = delete;
	BCryptAlgorithmProvider& operator=(const BCryptAlgorithmProvider&) = delete;
};

class BCryptHash : BCryptAlgorithmProvider
{
	int hash_length_;
public:
	int hash_length() const { return hash_length_; }

	BCryptHash(LPCWSTR algorithm, LPCWSTR implementation, ULONG flags);

	void ComputeHash(const void * buffer, size_t length, void * hash, size_t hash_length);
};

class BCryptSha256 : public BCryptHash
{
public:
	BCryptSha256() : BCryptHash(BCRYPT_SHA256_ALGORITHM, NULL, BCRYPT_HASH_REUSABLE_FLAG)
	{ }
};

// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "compute_descriptors.h"
#include <iostream>

void parallel::recursive_compute(const boost::filesystem::path & input_dir, int max_order, std::size_t queue_size, std::size_t max_thread, sqlite::database & db)
{
    using namespace std;
    using namespace boost::filesystem;
    using namespace logging;

    logger_t & logger = logger_main::get();

    TasksQueue all_voxel_paths{ queue_size };

    vector<thread> working_threads{ max_thread };

    atomic_bool is_stop{ false };

    using NodeType = tree::Node<std::string>;

    tree::PathTree<NodeType> tree{ "" };

    try
    {
        tree = tree::PathTree<NodeType>::build_tree_from_db(input_dir, db);
    }
    catch (const sqlite::sqlite_exception & exc)
    {
        BOOST_LOG_SEV(logger, severity_t::error) << "Terminate main thread." << endl << exc.what() << endl << exc.get_code() << endl << exc.get_sql() << endl;
        return;
    }

    for (size_t i{ 0 }; i < working_threads.size(); i++)
    {
        working_threads.at(i) = thread(compute_descriptor, ref(all_voxel_paths), max_order, ref(is_stop), ref(db));
    }

    auto iterator = recursive_directory_iterator(input_dir);

    // hex string
    string file_hash(picosha2::k_digest_size * 2, '\0');
    vector<unsigned char> hash_buffer(picosha2::k_digest_size, 0);

    for (const auto & entry : iterator)
    {
        if (is_stop)
        {
            break;
        }

        if (entry.status().type() == file_type::regular_file)
        {
            path local_file{ entry.path() };

            if (local_file.extension() == u8".binvox")
            {
                BOOST_LOG_SEV(logger, severity_t::info) << u8"Found " << local_file << endl;

                if (!::hash::compute_sha256(local_file, hash_buffer, file_hash))
                {
                    BOOST_LOG_SEV(logger, severity_t::warning) << u8"Cannot compute hash for " << local_file << endl;
                    is_stop = true;
                    continue;
                }

                shared_ptr<NodeType> node;

                bool is_new_node = tree.add_path(local_file, file_hash, node);

                path relative_path = relative(local_file, input_dir);

                bool need_recompute{ true };

                if (!is_new_node)
                {
                    if (file_hash != node->data())
                    {
                        BOOST_LOG_SEV(logger, severity_t::debug) << u8"File: " << local_file << " changed. Need to recompute." << endl;

                        stringstream delete_query;

                        delete_query << u8"DELETE FROM " << db::DbSchema::table_name() << " WHERE " << db::DbSchema::path_column() << " = ?" << relative_path.generic_string();
                        try
                        {
                            db << delete_query.str();
                        }
                        catch (const sqlite::sqlite_exception & exc)
                        {
                            BOOST_LOG_SEV(logger, severity_t::error) << exc.what() << endl << exc.get_code() << endl << exc.get_sql() << endl;
                            continue;
                        }
                    }
                    else
                    {
                        stringstream select_query;

                        long count{};

                        select_query << u8"SELECT count(*) FROM " << db::DbSchema::table_name()
                            << " WHERE " << db::DbSchema::path_column() << " = ? AND "
                            << db::DbSchema::max_order_column() << " = ?";

                        try
                        {
                            db << select_query.str()
                                << relative_path.generic_string()
                                << max_order
                                >> count;
                        }
                        catch (const sqlite::sqlite_exception & exc)
                        {
                            BOOST_LOG_SEV(logger, severity_t::error) << exc.what() << endl << exc.get_code() << endl << exc.get_sql() << endl;
                            continue;
                        }

                        if (count > 0)
                        {
                            need_recompute = false;
                        }

                        if (need_recompute)
                        {
                            BOOST_LOG_SEV(logger, severity_t::debug) << u8"Cannot find computed descriptor for: " << local_file << u8" when max_order = " << max_order << u8" Need recompute." << endl;
                        }
                    }
                }

                if (need_recompute)
                {
                    tuple<path, path, string> item = std::make_tuple(input_dir, relative_path, file_hash);

                    while (!all_voxel_paths.push(item) && !is_stop)
                    {
                        this_thread::sleep_for(500ms);
                    }
                }
                else
                {
                    BOOST_LOG_SEV(logger, severity_t::info) << u8"File: " << local_file << u8" with hash: " << file_hash << u8" and max_order = " << max_order << u8" already exists. Skip" << endl;
                }
            }
        }

        while (!(all_voxel_paths.empty() || is_stop))
        {
            this_thread::sleep_for(500ms);
        }
    }

    is_stop = true;

    for (auto & thread : working_threads)
    {
        thread.join();
    }

    BOOST_LOG_SEV(logger, severity_t::info) << u8"Completed" << endl;
}

void parallel::compute_descriptor(TasksQueue & queue, int max_order, std::atomic_bool & is_stop, sqlite::database & db)
{
    using namespace std;
    using namespace boost::filesystem;
    using namespace logging;

    using VoxelType = bool;
    using Container = vector<VoxelType>;
    using DescriptorType = double;

    Container binvox_voxels;
    Container canonical_order_voxels;
    size_t dim{};

    logger_t & logger = logger_main::get();

    tuple<path, path, string> path_to_voxel;

    using Row = sqldata::Row<DescriptorType>;

    const size_t rows_buffer_size{ 10 };

    sqldata::CollectionRows<DescriptorType> rows;

    while (true)
    {
        if (!queue.pop(path_to_voxel))
        {
            if (is_stop)
            {
                break;
            }

            std::this_thread::sleep_for(50ms);
        }
        else
        {
            path absolute_path = get<0>(path_to_voxel) / get<1>(path_to_voxel);

            BOOST_LOG_SEV(logger, severity_t::debug) << u8"Processing " << absolute_path << endl;

            if (!io::binvox::read_binvox(absolute_path, binvox_voxels, dim))
            {
                BOOST_LOG_SEV(logger, severity_t::warning) << u8"Cannot read binvox from " << absolute_path << endl;
            }
            else
            {
                canonical_order_voxels.resize(binvox_voxels.size());
                binvox::utils::convert_to_canonical_order(binvox_voxels.begin(), canonical_order_voxels.begin(), dim);

                // compute the zernike descriptors
                // This invoke changes voxels data
                ZernikeDescriptor<DescriptorType, Container::iterator> zd(canonical_order_voxels.begin(), dim, max_order);

                auto invs{ zd.get_invariants() };

                auto grid_ = vector<vector<vector<std::complex<DescriptorType>>>>(dim, vector<vector<std::complex<DescriptorType>>>(dim, vector<std::complex<DescriptorType>>(dim)));
                zd.Reconstruct(grid_);
                std::cout << "reconstruction is done \n";

                if (rows.size() < rows_buffer_size)
                {
                    rows.emplace_row(
                        get<1>(path_to_voxel).generic_string(),
                        get<2>(path_to_voxel),
                        invs,
                        max_order);
                }
                else
                {
                    try
                    {
                        db << rows;
                        BOOST_LOG_SEV(logger, severity_t::info) << u8"Save invariants to database." << endl;
                    }
                    catch (const sqlite::sqlite_exception & exc)
                    {
                        BOOST_LOG_SEV(logger, severity_t::warning) << u8"Cannot save invariants to database." << exc.what() << endl << exc.get_extended_code() << endl << exc.get_sql() << endl;
                        is_stop = true;
                        return;
                    }

                    rows.clear();
                    rows.emplace_row(
                        get<1>(path_to_voxel).generic_string(),
                        get<2>(path_to_voxel),
                        invs,
                        max_order);
                }
            }
        }
    }

    // Rest items
    if (!rows.empty())
    {
        try
        {
            db << rows;
            BOOST_LOG_SEV(logger, severity_t::info) << u8"Save invariants to database." << endl;
        }
        catch (const sqlite::sqlite_exception & exc)
        {
            BOOST_LOG_SEV(logger, severity_t::warning) << u8"Cannot save invariants to database." << exc.what() << endl << exc.get_extended_code() << endl << exc.get_sql() << endl;
            is_stop = true;
            return;
        }
    }
}

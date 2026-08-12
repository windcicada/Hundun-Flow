# Stage 3 final acceptance draft

status=ACCEPTED
accepted_governance_code_head=0cbd3d5bde4be63bc6346b4b32db771d87c59ea2
accepted_task11_head=66080e324089599711fdb26082af9b330bfdb5ce
candidate_parent=fe9065f8559e1367e8e112505bdd565f108d217f
candidate_tree=d50c1236f67bd2bdde58c94a125e530ae0f2ffea
task11_to_candidate_diff_sha256=6d244b21b16766a10039793b23eddc4af9e3ff97f70e82b8b257f1d536fbfac7
scientific_matrix_head=0cbd3d5bde4be63bc6346b4b32db771d87c59ea2
product_projection_source_tree=d50c1236f67bd2bdde58c94a125e530ae0f2ffea
accepted_product_head=22ed17b438ffbb121ccda97898580183bd0803f8
worktree=clean
background_processes=none
dco_after_activation=29/29
manifest_count=57/57
low_cost=16/16
sanitizer=10/10
governance=5/5
scientific=23/23
performance=3/3
capability_ledger_sha256=977e03373c4389d3fe0ddf0139fdd9d5b7d2c0202275b5a4fd6eb9763240b894

## Build roots

- /home/wyf/code_dev/.worktrees/hundun-flow-stage3-framework/build/stage3-final-debug/CMakeCache.txt cache_sha256=8c44af17600becf325e8e781fb354045c54ad130e877f6488d6c98e1653efdc8
- /home/wyf/code_dev/.worktrees/hundun-flow-stage3-framework/build/stage3-final-release/CMakeCache.txt cache_sha256=490ef7fc5d68be212017de0fcce92c876049791ddaa59ab7114513058471f89d
- /home/wyf/code_dev/.worktrees/hundun-flow-stage3-framework/build/stage3-final-tests-off/CMakeCache.txt cache_sha256=3db4354c7bcbecc52430213580d96632c827d7cc80a1993143014da2e362dbd2
- /home/wyf/code_dev/.worktrees/hundun-flow-stage3-framework/build/stage3-final-asan/CMakeCache.txt cache_sha256=501333e4226cf5d0114546ea0d9cde26e108809abb40c10c5dcc3bd8069b9cd3
- /home/wyf/code_dev/.worktrees/hundun-flow-stage3-framework/build/stage3-final-ubsan/CMakeCache.txt cache_sha256=e26e511858c73116513812dedcb29947bfd90934fc13ab5255631d8d3caebf71

## Toolchain

cmake=cmake version 3.31.12
cxx=/home/wyf/.local/bin/clang++
mpi=mpiexec (OpenRTE) 2.1.1

## Commands

- five-root V0 configure/build matrix: exit 0
- stage3_acceptance --group low-cost: 16/16 exit 0
- stage3_acceptance --group sanitizer: 10/10 exit 0
- stage3_acceptance --group governance: 5/5 exit 0
- stage3_acceptance --group scientific: 23/23 exit 0
- stage3_acceptance --group performance: 3/3 exit 0
- manifest identity audit: 57/57, 0 errors

## Terminal manifests

| row | group | exit | manifest_sha256 | log_sha256 | binary_sha256 |
| --- | --- | ---: | --- | --- | --- |
| asan-checkpoint-density-profiles-fast | sanitizer | 0 | 6a511912a84f2abc84bfaec72202c1757ced30bfd3999b1ee082bfcc77bc61f9 | 9690ede6a1c50382e4f86f1646eacfaa98be2e8b9b3a6d092ff5577d9e46810a | df355c31910c5817d3a0793c122705cb718e0c884af742262372f8ea94905c5c |
| asan-evidence-manifest-unit | sanitizer | 0 | 0da450ecd54c04e1ea128c93bcbb44ca8012f10941e9067d755cfb205be9cd8a | f3cbf618ea495cf34d105d23b2088509589205c1861bd5b8c46a91334627258a | f0feadbebe8fde538b8b9d5436d78860c50c4c5178968371fb41ebda784d7798 |
| asan-ideal-wale-composition-fast | sanitizer | 0 | c7cb7f15cc0437342748855f49aa36e412ffb77eed989a092cdf10eff2e36459 | 87d069b05a6bee3bb4f184712d5385455052776768350632f7ca2ba7157e7b89 | 268df846cf3e4d22397622b31d6137e22b70ef2f132dcfee067df42d6d1faaee |
| asan-immersed-wale-constant-fast | sanitizer | 0 | 523704be1e024596b1d5b68fd494a8bd637ca71f0e3ea9543c6f1bd8a46e8811 | 143432f9081f7809286caae8583178c41f6c4b77169763185cebafc3e92b92c3 | e09b690207f545d5baaaf6e297004d31df9a2fd6d360290a54e9f8118fac22f1 |
| asan-material-wale-composition-fast | sanitizer | 0 | e1048b332668f4b02766ab1ce6887272c9003b8a512f80cc345ac472f5ab3d54 | e8194195fc7ecebc142683b671bcb7191128eafdf1abc432f519c9800ae9f2c0 | f6f28905e1bdf4362ae5019fde663c6bd3f831c8404da8633c39b6e71c16c841 |
| checkpoint-continuation-n12-r1 | scientific | 0 | d6c11f81c71bab4ef79387b661f04576b47e7ebd9d584f02e7a1dc11b1e8d70a | 3402c368b801219ebbe16f6a59e3a47182fdcc28b214a549b51728f4f1a351d3 | c4b45f693a3ea93994c274903eb123cd36f14af5e65bdc3397dc923a25b93c34 |
| checkpoint-continuation-n12-r2 | scientific | 0 | c5869370e7bee606013a44e2c6d8b77cdc7c20f0ad39cb39a183c641cd304949 | 783538cc5203224085174390b7528b13edb5121c5cb57fe1d388c29d77623678 | c4b45f693a3ea93994c274903eb123cd36f14af5e65bdc3397dc923a25b93c34 |
| checkpoint-continuation-n12-r4 | scientific | 0 | 17445baaf61621bb0c8e90dda20a44941dfc4cecbd8788b46c7bcd79a822ef87 | cba90b4e0fb7bd288b98b7154928a8795815e1be93e4b5f72a95e1ded6fe59d5 | c4b45f693a3ea93994c274903eb123cd36f14af5e65bdc3397dc923a25b93c34 |
| constant-ibm-wale-n24-r1 | scientific | 0 | 3c4deaa7dec0d822eb2b9e612b53c7b5e88d35546e69722178d2b336348ffdb7 | 6866769afd5a6fd5a502e66c341c8bac85bbff6f31d7f97f1a9d5e6e607c838d | 91fb4ac5fcf97b6271121809d8cf68243490cad4d4f0b3ab982ba983f54c963b |
| constant-ibm-wale-n24-r2 | scientific | 0 | 3bf89b43896909c34370abf2b84b298164dc8678cd2c404e45be65bd00267d4a | 8a2540b23184e7c5e0bc09713e4122144fbc3568bbf2464f03c9212c296b0274 | 91fb4ac5fcf97b6271121809d8cf68243490cad4d4f0b3ab982ba983f54c963b |
| constant-ibm-wale-n24-r4 | scientific | 0 | e270f603e8f1c7d4c584ce4e6778a83ed100ca4329e634d98753554c3990ec07 | 5b6c8feb8dbf2f1fc2450fd6f97fc19f78732731da568689166bd53670702e5c | 91fb4ac5fcf97b6271121809d8cf68243490cad4d4f0b3ab982ba983f54c963b |
| constant-ibm-wale-n48-r1 | scientific | 0 | 8e5dcea6ca6c83522de2486d7f0d5c463c2d573a9abd0c0ae334c3fed068aaef | eb08325ccff01b931505fa32b6165505cf5cfc77941b571e1e1daf35cca59a59 | 91fb4ac5fcf97b6271121809d8cf68243490cad4d4f0b3ab982ba983f54c963b |
| diagnostics-fast-r1 | low-cost | 0 | a0bf26e7ddf71a5f508a1088eade076d681f6967c94e46a470bd6e04b6bdbcea | 5c96dc7156bff4bc5d853ff6ebe231c127bcd794a1fafe7e7596403db16a3b43 | 151dd63f0c42722e266ce8dcd3fc28bcb0c919f5c4cea9fa53282100647cc6ed |
| diagnostics-fast-r2 | low-cost | 0 | fd46cb1c1e2876dd9b2af14da553e260ff97d5fb6454cd918c483290aa9b76cb | fbe4cb849b45e682f1eec36117f3858b2a6f9ed149b259ac5d6332030f57d659 | 151dd63f0c42722e266ce8dcd3fc28bcb0c919f5c4cea9fa53282100647cc6ed |
| diagnostics-fast-r4 | low-cost | 0 | 759e8e467facb0835614c5817e29b9aeeea7a64c414666f242602984e75719fd | f6cd35e8d24a53c8d699430f9038aa5bd31e1b0618298150f2e122532be20766 | a3a0ab85881e9f8cad225f3e2e7bcd5a2ef19911dc84e33abd1e80d2335986cc |
| governance-ledger | governance | 0 | a6c982029046c5ee0806dd552c7403e0c4d37c34d9bee648a78e406c34d7bb04 | dbea4fcefdac59dd01419e0db0ff5bd5226ef1dbb46b713d0dde7481f25b7e60 | 8c44af17600becf325e8e781fb354045c54ad130e877f6488d6c98e1653efdc8 |
| governance-projection | governance | 0 | 387d84ec6b9636d976c32a694f5b74d0ae0829a88ba11b0285c25605b82649ab | fb804d1465cbaa4dfe62627b7eb8cbc380cdf3fbe612b5a6a3587f5a3d377eae | 8c44af17600becf325e8e781fb354045c54ad130e877f6488d6c98e1653efdc8 |
| governance-provenance | governance | 0 | 39cdd5591ae9c9ff65ca395864f41b04035ed8b88d2a263c0dabc96ee4665edc | a284097b267d6d2b41dd40068879140d7f96b8b04f41389480b26602f9713643 | 8c44af17600becf325e8e781fb354045c54ad130e877f6488d6c98e1653efdc8 |
| governance-source-policy | governance | 0 | 05ade0a6eb734542bddb81bec7965f5b9c22c3ca034c33c7f8b366e3f9cdf65e | 888724ab1f818e3ab29e4fba6a98dc2ed2345bfbf3edaa5529bc662ea50044f3 | 8c44af17600becf325e8e781fb354045c54ad130e877f6488d6c98e1653efdc8 |
| governance-tests-off | governance | 0 | 01646d98e8c1a3c32fa5da28ed57281e11f08c29e1620a8267752afae8b89e1f | 643a92701963062cacf347a5a1f0502b7e506e72bfed47cc2398faf8319dc005 | 996269ee9f6eea3ebfd5c0d5ade5132153dbb6124c1d0c3851dcfe4fad0eb9ab |
| ideal-ibm-wale-n12-r1 | scientific | 0 | 876327673ab7e9ca5f0ac11346dab5385ba0e14943f7964990228a403b8b0db4 | 79f56f2112c15641613708171a973a8136d7773bd6e7b0dfd89a484bb651664a | 047f7972f66aa4343add6717a13fcead751457ed28f41bbdd77ee94c1dfc30d1 |
| ideal-ibm-wale-n12-r2 | scientific | 0 | c11be57995a5510b52acf594b207efcb5da106208ede0658d882a9b9425109fb | e7ffa4605b377848375bed0fc01be7ba337ad749828a9a9931d28428b01da32a | 047f7972f66aa4343add6717a13fcead751457ed28f41bbdd77ee94c1dfc30d1 |
| ideal-ibm-wale-n12-r4 | scientific | 0 | 478df0e2648d3610a5c0959ec09aec1da2ed33e97a951b2daee08adfd1fdeb52 | 938d13d69012a1fe1a80a95b3a05014a3e853d4352945ce1d9b444282034ec38 | 047f7972f66aa4343add6717a13fcead751457ed28f41bbdd77ee94c1dfc30d1 |
| ideal-ibm-wale-n24-r1 | scientific | 0 | 226a08ee041a52351f3289d41e996e3c58dcdcf24d888b852eef74864e9dede2 | 34d963a99308bdb8812745e131b6864d78f73da5f83596347a7ed6150abba79c | 047f7972f66aa4343add6717a13fcead751457ed28f41bbdd77ee94c1dfc30d1 |
| ideal-ibm-wale-n24-r2 | scientific | 0 | 539178db14eda0618526cd36924d5fb17b5b5a93f0445e16d2e60c68098baa04 | 1a540d64e413fac43ebc192930c431461069f6c101bccab85f5a64d31754ba38 | 047f7972f66aa4343add6717a13fcead751457ed28f41bbdd77ee94c1dfc30d1 |
| ideal-ibm-wale-n24-r4 | scientific | 0 | 8bd99d04bc3105aed5e0378cff141e999463406a78f39238f782b5e66ed0853f | 8d66543cf0a751570ce44f6247c5eac268688bd3a101cb67edf26aec376ddf59 | 047f7972f66aa4343add6717a13fcead751457ed28f41bbdd77ee94c1dfc30d1 |
| material-ibm-wale-n12-r1 | scientific | 0 | ac2dc5ba6e3a62c4dd0ecdb5cd75e425577181d99b520ee336a48e4d8e0820ec | 69df3dffb0cab73a881c4971af6eda2c09a27080f836f719ebb207db44f53ec6 | 175d9aa2309356b3a76d3eeed7252e0108a240e5c8a70b8b5fd37693208e2736 |
| material-ibm-wale-n12-r2 | scientific | 0 | 3f0d605e4f3b27989ae70ef75d11097a45fa7a96e4f75aec5abf47a0073ca0f2 | efe5664cade5253999701a7a95fe434a812aefc53662df17651c0155a3926e5d | 175d9aa2309356b3a76d3eeed7252e0108a240e5c8a70b8b5fd37693208e2736 |
| material-ibm-wale-n12-r4 | scientific | 0 | 1ed6f25e1bf20d1b23708a1e5801249e24a8fb37c019d71a5d831cac641697f7 | 7efc602e9b579609534905708bd4a6e84aaee333a05fad9920ab338b87dd968a | 175d9aa2309356b3a76d3eeed7252e0108a240e5c8a70b8b5fd37693208e2736 |
| material-ibm-wale-n24-r1 | scientific | 0 | 80070a287296c136f15d9033a84d25bbf5dca8f59a209890a42f15232df33144 | f761b5b920cf6a4195fabfa667cdd40376669d19ad079d7d251806bfed284c93 | 175d9aa2309356b3a76d3eeed7252e0108a240e5c8a70b8b5fd37693208e2736 |
| material-ibm-wale-n24-r2 | scientific | 0 | 9863dea7571da40983cd81ddcc18e4bc8272452b1fbe3ba74af8fcddffdc1953 | f77f1a2669fe2f027dc2f8565b43d16f28b3eb391a7c13e7121fca59a2abdf12 | 175d9aa2309356b3a76d3eeed7252e0108a240e5c8a70b8b5fd37693208e2736 |
| material-ibm-wale-n24-r4 | scientific | 0 | 5a966bd3b505b3589fa95743422506fb31700b67ca9a6dcca21bde7ce90e0315 | a91917502cf4f77f96d4d1dab39e2806a63c6a8b2f088679223c2d271c4a516d | 175d9aa2309356b3a76d3eeed7252e0108a240e5c8a70b8b5fd37693208e2736 |
| profiles-fast-r1 | low-cost | 0 | d4dec1c0288080317a3e93a957d868b77f9cda4fe60bb42bbcc1a2f9dad6f649 | 9c52290577dc0fa89e51f6c884380fad5b0736b107d7ef82d508e5d0d44ffa99 | 8c44af17600becf325e8e781fb354045c54ad130e877f6488d6c98e1653efdc8 |
| profiles-fast-r2 | low-cost | 0 | 22799df6fc8657712d67cebae8a54f7bf88ec9bf2ded0e2fb61b742370373033 | 25d0b20801266ee9916ff68c3f625cd4a1dfa1a0abfa73c550d5b1f318f87836 | 8c44af17600becf325e8e781fb354045c54ad130e877f6488d6c98e1653efdc8 |
| profiles-fast-r4 | low-cost | 0 | 393eae4065e38c9b9544ece605ff3b10f72d4b714250a4fe4fe4bf4554f690f8 | 18b0d78f9061edc8afdfa9a36337dfe2243f7ba30f77d2114813ab511647c76f | 8c44af17600becf325e8e781fb354045c54ad130e877f6488d6c98e1653efdc8 |
| registration-contract | low-cost | 0 | 0dec34300d12ef8001aa205a30cdfc98cc651dbd896cfd340a9bb5c009873803 | aeb5b42f74f5aca50c314e101032c64030ee83e814d5aa4f9b0eba616b6a1617 | 8c44af17600becf325e8e781fb354045c54ad130e877f6488d6c98e1653efdc8 |
| restart-fast-r1 | low-cost | 0 | 36b9d12054728c3a9fb9553b85ccf6275d152ddaa16c38d3283372dc7e99a151 | d3c72f15af3bd09c38c49d8e9534e0298a368ac9e64a04806c16ecdda551f650 | 8c44af17600becf325e8e781fb354045c54ad130e877f6488d6c98e1653efdc8 |
| restart-fast-r2 | low-cost | 0 | 13f2df46fb8a0df2b74d62eded23b3e6a80b5fa939ae8dbfa90341c893f1a908 | 63a6bedc7c15dfb179f996b4006112407d3bf1d445942d34e4336d7d55ce2392 | 8c44af17600becf325e8e781fb354045c54ad130e877f6488d6c98e1653efdc8 |
| restart-fast-r4 | low-cost | 0 | ad598cfa6815ab136c97f39d7f87dc68160f45854e510b440fdd0d01619c752f | fccfde18974141c358886eac872d8015a5cebda63ab6a3a9b0f276a4ad1559f5 | 5cf6ae8c7ee772fbd5912005ef07cb15964ebbe12d4e7ded862a333e0bb5afce |
| science-row-contract | low-cost | 0 | aef640099476102380ff0acdef8d03ce8e5b8e598bdc2b6537836084e054b3a5 | 9cd607e73351e08c124e857ca2834df98b7ab98e7eb8994b71a5cb2e6d655653 | 1157b1be86d49ee20985e04c188763544d1640c2e9df36808cb40d1d89df0ed7 |
| stage1-low-cost-whitelist | low-cost | 0 | b7cb37978fe9d51d39d95f946f2881f20f706c5193b652716cd4512d17ab309e | 6235392784743e5aaa7a455df0cd38f15c13703342426ceff2a4d7dc8cada8b2 | 17063ed824690e3f171d1609e7297bc7b98319a0d8717307c55947835b02d7d6 |
| stage2-core-whitelist | low-cost | 0 | d93e242202de6505ff004530b54ffde9f038f51dd072bf3af4b55dccd54f7ea4 | 2800d996277c1f0e80186b0fc310ec49928cdd546c99bd2db23dbe5b598bb1bc | 8c44af17600becf325e8e781fb354045c54ad130e877f6488d6c98e1653efdc8 |
| stage3-exact-counters-fast-r1 | low-cost | 0 | a2f69b97e85228ff7937890b42ccca1269cb47902cdfb4f9390fcd8bc692c397 | 6e79b83dd8b89b761d790d696b28c43239d713dfa9870df3d23335c30849de24 | aa10deaf91df95664b6ba8158c70c192c10668a71305e62fdbbccf9a33408ddd |
| stage3-exact-counters-fast-r2 | low-cost | 0 | 5e40a30ea6b0e1f62cc5580fcf0124532a88a04c866922f37ba056eb95bc1064 | 813836e845e858819626ca3094bfe50c45b21b077ed4b4bf0ecfdf4c2064a13c | aa10deaf91df95664b6ba8158c70c192c10668a71305e62fdbbccf9a33408ddd |
| stage3-exact-counters-n24-r1 | performance | 0 | b3b6d391d3b93c463a5a6562acff816b3f8d2467fe1940a0d57a6b39039a7f1f | fa1ae0b187cb7d0952d99635ba244a341db3a1b7033d7efe3a6b1287eee04e93 | 5213d24211c457655d45f8d5d304103fdbdb7cbdbf8e9f2f891a326d6e071a27 |
| stage3-exact-counters-n24-r2 | performance | 0 | 7e5c3b1b0446436344c731c6ed97a6d99a17871680af6a71b2d18fccb517d706 | 438aa6314db7c12d7d66f9c8e50803d60feca2b4d2b4fb14433b778e68c7cf52 | 5213d24211c457655d45f8d5d304103fdbdb7cbdbf8e9f2f891a326d6e071a27 |
| stage3-exact-counters-n24-r4 | performance | 0 | a608e4e34f08a82b4803043dd48c9a7089904af0fa35f4270c130507b47e141d | 349eff137d2c57a966a633bdaab9fabea2a2cb343f6999ed4f908dc100047545 | 5213d24211c457655d45f8d5d304103fdbdb7cbdbf8e9f2f891a326d6e071a27 |
| task11-authority-current-tree | low-cost | 0 | bda81d401ec745da0678ca97f30849d079eaca5b6ec8f7c4331b4b76be844318 | 8dbc864e9dc6820d2955dbd0413a60950a0c486b0ead17023e3af3550d80f843 | 8ebedeecbaf139f7d948bf29afd53579d48d7a8f7c867ec01deed60420684977 |
| ubsan-checkpoint-density-profiles-fast | sanitizer | 0 | 73bcd7007d867cbc664cb6826b893912f5f2b4649a665a1fcc1377a141c39e40 | 4e3442d0720d31467038815a57b50c5be25aebf9e582768ab82e46c876427554 | 59bc84ce371aded30dbbb4ade2483bb5bb03c5448a82a51bdb3f97cb893f2e53 |
| ubsan-evidence-manifest-unit | sanitizer | 0 | 10e1450277dfbe70e5094639c0a943070b422057911b22bd1d55e0633d9cc17d | deb4f70023ecfe70ccde5f01086e48f2a01537cfcb9271b21f60a134ee12d4e2 | b9c86a1353eacb9e0d8056aa07849aa9bd7fd21c16cd1533f4983f9d4cb9b670 |
| ubsan-ideal-wale-composition-fast | sanitizer | 0 | e94ef9b13dc1114ce251aa81583dda7fe5df208df6ab94b152bd484adbbe8cdd | c18e06f8059a71caac4af35600620ed86e3d7eb4eb92f82d89081fd627937a3e | a34736b6fc8da1996f11ab5a695e59def10e0a2738eb3ce41df36851b22466b4 |
| ubsan-immersed-wale-constant-fast | sanitizer | 0 | a4f99a6b64f61e2e28ecf18e87ca1e0df900205de221933436abe1f51967cae3 | 492523e93bcce65f06957ea5c8ab097c6a6716220a02098f261ddd4a4d8949b4 | 9ca5fab5e8b1d1852f8b738443ce426899f6625ba97f482415e703cd2df0d2b0 |
| ubsan-material-wale-composition-fast | sanitizer | 0 | 3ca63a37ca1bf3139033d0161276ea19166d9593efbf64cd7bbd76e3a3e9de3a | 0fc7009fb6c7e2c29871c0d0f50424a48a953995cbd2fa326958816b7b24127b | e4246c8907de1b3b51991b64628e635575f4f032da5061e94c82bc5cb6d02cf5 |
| wale-channel-n48-r1 | scientific | 0 | 9daa6c37f7eee810693c469c5997abe3281c8b70e7dfaae9cb5a224018a4238c | fd30bc54cb3819beec9c30f71879f0a00cb027c76861da29c6a7f36e2dbdee31 | cae5a11d26ed56f251e52abf749b2f7b105b731d06a52b1c49e8918bfdda929b |
| wale-tgv-convergence-r1 | scientific | 0 | 60ee31b53b127e8c262005ffb8326c78be1a540cfeb85a88f3895e081116b578 | 466b26b5472b03847b119a4bc21f83b2b901b375cbe200b54237f3ed8b8c762f | c5d6e13950d8195010a85049a6660e6050f347c220b142532bc7fca3c7b8d175 |
| wale-tgv-n24-r2 | scientific | 0 | f3ed05cdd4823f1cba07d2c4fd7e6e7cf32c49a820d50994c189b24d9f4103cd | 17ac6c6807ad14347faf5777b35b8810e73cdac0be787a84774e24df0f64e926 | c5d6e13950d8195010a85049a6660e6050f347c220b142532bc7fca3c7b8d175 |
| wale-tgv-n24-r4 | scientific | 0 | ca41703ea259f7843cea727a02283b2eec2906a13c2c777110d659e5dd4faa5c | 669c081e547504b0015586897008c1c6d524a2cf1306399c837b4d1548dba666 | c5d6e13950d8195010a85049a6660e6050f347c220b142532bc7fca3c7b8d175 |

## Product projection

- base_product_head=ae3d08bbb220d1d3b28ec070d1cba9c33fb85877
- accepted_product_head=22ed17b438ffbb121ccda97898580183bd0803f8
- accepted_product_tree=7fb9ce848238eeab5dc1ad0908092d8d115851b4
- product_projection_manifest_sha256=224a3cdbb6fb104ad103256eb0de28732ff5e9dfbe47f8b120460bac2ea25f8c
- product_projection_paths=272/272
- product_build_binary_sha256=7aaaf90e34662c19e300f116703943a0c4559feea94054b2d92c9a2a83152bf9
- product_installed_binary_sha256=742a1fea5cc562a790e9bd3d6809498e50869f514bfbf6a33c0bd0fd9ee9fa1e
- product_root_fast_forward=PASS
- product_remotes=none
- product_worktree=clean

## V2 command record

- manifest-bound projector: 272 paths, exit 0
- pre-commit product path/blob/text/license/NOTICE scan: 272/272, exit 0
- default GCC 7.5 diagnostic build: exit 2, rejected because the environment compiler lacks the required C++17 filesystem header
- explicit Clang without libc++ diagnostic build: exit 2, rejected because the incomplete toolchain omitted the authenticated standard-library setting
- Clang+libc++ hundun-only build: exit 0
- initial install after hundun-only build: exit 1, rejected because the full install set required the not-yet-built SDK archive
- full Clang+libc++ tests-off build/install: exit 0
- built and installed hundun --version: HUNDUN-FLOW 0.2.0, exit 0
- built and installed template --validate: VALID, exit 0
- built and installed template --print-resolved: identical, exit 0
- installed one-step 1-rank smoke: FINISHED step=1 time_s=0.001, exit 0
- staged product audit: 272/272 with 80 changed paths, exit 0
- signed product commit: 22ed17b438ffbb121ccda97898580183bd0803f8
- post-commit current-tree and complete-history scan: 2 commits, 272 paths, zero remotes, exit 0
- product-root fast-forward: ae3d08bbb220d1d3b28ec070d1cba9c33fb85877 -> 22ed17b438ffbb121ccda97898580183bd0803f8, exit 0

## Dispositions

- 96-cubed: permanently excluded and not run.
- AMR, GPU, moving bodies, rank-changing Restart: deferred/out of Stage 3.
- Private source/data access: none.
- Research-process interference: none.
- Push/publication: none.
- Product commit and root fast-forward: accepted at 22ed17b438ffbb121ccda97898580183bd0803f8.
- Tracked governance seal: this three-file report transaction; its commit ID is reported after commit.

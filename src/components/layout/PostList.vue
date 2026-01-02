<!-- PostList.vue -->
<template>
  <div>
    <div v-for="post in displayedPosts" :key="post.id" class="col-12 mb-3">
      <div class="card shadow-sm border-0" :style="{ backgroundColor: 'var(--app-card-bg)' }">
        <div class="card-body p-4 d-flex flex-column flex-md-row gap-4">
          <div class="flex-grow-1">
            <small class="meta-text mb-2 d-block" :style="{ color: 'var(--app-text-muted)' }">
              {{ formatDate(post.date) }}
            </small>

            <h5 class="post-title fw-bold mb-2" :style="{ color: 'var(--app-text)' }">
              <router-link :to="getArticlePath(post)" class="text-decoration-none"
                :style="{ color: 'var(--app-text)' }">
                {{ post.title }}
              </router-link>
            </h5>

            <div class="mb-2">
              <small class="meta-text" :style="{ color: 'var(--app-text-muted)' }">
                {{ formatCategory(post.category) }}
              </small>
            </div>

            <p class="mb-3 desc-text">
              {{ post.preview }}
            </p>


            <div class="d-flex flex-wrap gap-2">
              <span v-for="tag in post.tags" :key="tag"
                class="badge tag-badge fw-normal py-1 px-2 rounded-3 cursor-pointer" @click="goTag(tag)">
                # {{ tag }}
              </span>
            </div>
          </div>
        </div>
      </div>
    </div>

    <div class="row mt-4 pb-4" v-if="totalPages > 1">
      <div class="col-12">
        <nav :aria-label="paginationLabel">
          <ul class="pagination justify-content-between align-items-center mb-0">
            <li class="page-item" :class="{ disabled: currentPage === 1 }">
              <button class="page-link d-flex align-items-center border-0 bg-transparent px-3 py-2"
                :style="{ color: 'var(--app-text-muted)' }" @click="prevPage" :disabled="currentPage === 1"
                :aria-label="t('prevPage')">
                <span class="me-1">&lt;</span>
                <span>{{ t('prevPage') }}</span>
              </button>
            </li>

            <li class="page-item">
              <div class="d-flex gap-2">
                <button v-if="showFirstPage" class="page-link d-inline-flex align-items-center justify-content-center"
                  :class="{ 'current-page': currentPage === 1 }" @click="goToPage(1)">
                  1
                </button>

                <span class="d-flex align-items-center" v-if="showFirstEllipsis">...</span>

                <button v-for="page in middlePages" :key="page"
                  class="page-link d-inline-flex align-items-center justify-content-center"
                  :class="{ 'current-page': currentPage === page }" @click="goToPage(page)">
                  {{ page }}
                </button>

                <span class="d-flex align-items-center" v-if="showLastEllipsis">...</span>

                <button v-if="showLastPage && totalPages > 1"
                  class="page-link d-inline-flex align-items-center justify-content-center"
                  :class="{ 'current-page': currentPage === totalPages }" @click="goToPage(totalPages)">
                  {{ totalPages }}
                </button>
              </div>
            </li>

            <li class="page-item" :class="{ disabled: currentPage === totalPages }">
              <button class="page-link d-flex align-items-center border-0 bg-transparent px-3 py-2"
                :style="{ color: 'var(--app-text-muted)' }" @click="nextPage" :disabled="currentPage === totalPages"
                :aria-label="t('nextPage')">
                <span>{{ t('nextPage') }}</span>
                <span class="ms-1">&gt;</span>
              </button>
            </li>
          </ul>
        </nav>
      </div>
    </div>
  </div>
</template>

<script>
/*
  PostList
  - 列表与分页组件，采用 props docs 与 perPage 控制
*/
import { useI18n } from 'vue-i18n'
import { loadNotes, loadCategories } from '@/utils/contentLoader'

export default {
  name: 'PostList',
  setup() {
    const { t, locale } = useI18n()
    return { t, locale }
  },
  props: {
    docs: {
      type: Array,
      required: true,
      default: () => []
    },
    perPage: {
      type: Number,
      default: 6
    }
  },
  data() {
    return {
      currentPage: 1,
      maxVisiblePages: 5,
      notesFlat: [],
      categoriesData: []
    }
  },
  computed: {
    paginationLabel() {
      return this.t('paginationLabel')
    },
    prevPageText() {
      return this.t('prevPage')
    },
    nextPageText() {
      return this.t('nextPage')
    },
    tagsText() {
      return this.t('tags')
    },
    totalPages() {
      return Math.max(1, Math.ceil(this.docs.length / this.perPage))
    },
    displayedPosts() {
      const start = (this.currentPage - 1) * this.perPage
      const end = start + this.perPage
      return this.docs.slice(start, end)
    },
    allVisiblePages() {
      const pages = []
      const total = this.totalPages
      const current = this.currentPage
      const maxShow = this.maxVisiblePages

      if (total <= maxShow) {
        for (let i = 1; i <= total; i++) {
          pages.push(i)
        }
        return pages
      }

      pages.push(current)

      let i = 1
      while (pages.length < maxShow) {
        if (current - i >= 1 && !pages.includes(current - i)) {
          pages.push(current - i)
        }
        if (pages.length >= maxShow) break

        if (current + i <= total && !pages.includes(current + i)) {
          pages.push(current + i)
        }
        i++
      }

      if (!pages.includes(1)) {
        pages.pop()
        pages.push(1)
      }
      if (pages.length < maxShow && !pages.includes(total)) {
        pages.push(total)
      } else if (pages.length >= maxShow && !pages.includes(total)) {
        pages.pop()
        pages.push(total)
      }

      return pages.sort((a, b) => a - b)
    },
    middlePages() {
      return this.allVisiblePages.filter(page =>
        page !== 1 && page !== this.totalPages
      )
    },
    showFirstPage() {
      return this.allVisiblePages.includes(1)
    },
    showLastPage() {
      return this.allVisiblePages.includes(this.totalPages)
    },
    showFirstEllipsis() {
      return this.totalPages > this.maxVisiblePages && this.allVisiblePages[1] > 2
    },
    showLastEllipsis() {
      return this.totalPages > this.maxVisiblePages &&
        this.allVisiblePages[this.allVisiblePages.length - 2] < this.totalPages - 1
    },
    categoryTitleMap() {
      const map = {}
      try {
        (this.categoriesData || []).forEach(section => {
          (section.items || []).forEach(item => {
            (item.categories || []).forEach(cat => {
              if (cat && cat.key && cat.title) map[cat.key] = cat.title
            })
          })
        })
      } catch (e) {
        console.error(e);
      }
      return map
    }
  },

  watch: {
    docs() {
      this.currentPage = 1
    },
    '$route.query.page'(newVal) {
      const p = parseInt(newVal);
      const page = Number.isFinite(p) && p >= 1 ? Math.min(p, this.totalPages) : 1;
      if (page !== this.currentPage) {
        this.currentPage = page;
        this.$nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }));
      }
    },
    locale() {
      this.loadData()
    }
  },
  async created() {
    await this.loadData()
  },
  mounted() {
    const p = parseInt(this.$route.query.page);
    this.currentPage = Number.isFinite(p) && p >= 1 ? Math.min(p, this.totalPages) : 1;
    this.handleResize();
    window.addEventListener('resize', this.handleResize);
  },
  methods: {
    formatDate(dateString) {
      const locale = this.locale === 'zh-CN' ? 'zh-CN' : 'en-US'
      const options = { year: 'numeric', month: 'long', day: 'numeric' }
      return new Date(dateString).toLocaleDateString(locale, options)
    },
    getArticlePath(post) {
      let articlePath = '';

      const found = Array.isArray(this.notesFlat) ? this.notesFlat.find(item => item.title === post.title) : null
      if (found && found.relativePath) {
        articlePath = `notes/${found.relativePath}.md`
      } else {
        const locale = this.locale;
        const isEnglish = locale === 'en-US';

        const basePath = post.category?.[1] || 'notes';
        const fileName = post.title.toLowerCase().replace(/[^a-z0-9]/g, '-');
        articlePath = `${basePath}/${fileName}.md`;

        if (isEnglish) {
          articlePath = articlePath.replace('.md', '-en.md');
        }
      }

      return `/article/${articlePath.replace(/\.md$/, '')}`
    },
    goToPage(page) {
      if (page >= 1 && page <= this.totalPages) {
        this.currentPage = page;
        const q = { ...this.$route.query, page: String(page) };
        this.$router.push({ path: this.$route.path, query: q }).catch(() => { });
        this.$nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }));
      }
    },
    prevPage() {
      if (this.currentPage > 1) {
        const page = this.currentPage - 1;
        this.currentPage = page;
        const q = { ...this.$route.query, page: String(page) };
        this.$router.push({ path: this.$route.path, query: q }).catch(() => { });
        this.$nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }));
      }
    },
    nextPage() {
      if (this.currentPage < this.totalPages) {
        const page = this.currentPage + 1;
        this.currentPage = page;
        const q = { ...this.$route.query, page: String(page) };
        this.$router.push({ path: this.$route.path, query: q }).catch(() => { });
        this.$nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }));
      }
    },
    goTag(tag) {
      if (!tag) return;
      const q = { ...this.$route.query, tag: tag, page: '1' };
      this.$router.push({ path: '/', query: q }).catch(() => { });
      this.$nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }));
    },
    clearTag() {
      const q = { ...this.$route.query };
      delete q.tag;
      q.page = '1';
      this.$router.push({ path: '/', query: q }).catch(() => { });
      this.$nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }));
    },
    handleResize() {
      this.maxVisiblePages = window.innerWidth < 480 ? 3 : 5
    },
    formatCategory(catArr) {
      if (!Array.isArray(catArr) || catArr.length === 0) return ''
      const top = catArr[0] || ''
      const subKey = catArr[1]
      if (!subKey) return top
      const subTitle = this.categoryTitleMap[subKey] || subKey
      return `${top} / ${subTitle}`
    },
    async loadData() {
      try {
        const [notesData, categoriesData] = await Promise.all([
          loadNotes(),
          loadCategories()
        ]);
        this.notesFlat = notesData || [];
        this.categoriesData = categoriesData || [];
      } catch (error) {
        console.error('Failed to load data:', error);
        this.notesFlat = [];
        this.categoriesData = [];
      }
    }
  },
  beforeUnmount() {
    window.removeEventListener('resize', this.handleResize)
  }
}
</script>

<style scoped>
.post-title {
  font-size: 1.4rem;
  font-weight: 700;
}

.desc-text {
  font-size: 1rem;
  line-height: 1.8;
  color: var(--app-text);
}

.meta-text {
  font-size: 0.95rem;
}

.badge {
  font-size: 0.85rem;
  font-weight: 500;
}

.tag-badge {
  color: var(--app-tag-text) !important;
  background-color: var(--app-tag-bg) !important;
}

.page-link {
  transition: all 0.2s ease;
  color: var(--app-pagination-link-text);
  background-color: var(--app-pagination-link-bg);
  border: none;
  outline: none;
  font-weight: 500;
  min-width: 36px;
  height: 36px;
  border-radius: 8px;
  display: flex;
  align-items: center;
  justify-content: center;
}

.page-link:hover {
  background-color: var(--app-primary-bg-subtle);
  color: var(--app-primary);
}

.current-page {
  background-color: var(--app-primary) !important;
  color: var(--app-pagination-current-text) !important;
  font-weight: 600 !important;
}

.page-link:focus {
  box-shadow: none;
  outline: none;
}
</style>

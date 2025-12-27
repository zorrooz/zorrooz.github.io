<!-- Resource.vue -->
<template>
  <div class="container py-4 px-3 view-container resource-view">
    <div class="row justify-content-center">
      <div class="col-lg-10 col-xl-8 typography-body">
        <div class="text-center mb-4">
          <h1 class="article-title mb-3">{{ pageTitle }}</h1>
          <p style="color: var(--app-text-secondary); margin-top: 0.75rem; margin-bottom: 0;">
            {{ pageSubtitle }}
          </p>
        </div>

        <div class="d-flex flex-column" style="gap: 2rem;">
          <section v-for="category in resources" :key="category.title">
            <h2 class="h4 fw-semibold pb-2 mb-3 heading-underline" style="color: var(--app-text-emphasis);">
              {{ category.title }}
            </h2>

            <div class="ms-3">
              <div v-for="sub in category.children" :key="sub.title" class="mb-4">
                <h3 class="h5 fw-semibold mb-3" style="color: var(--app-text-muted);">
                  {{ sub.title }}
                </h3>

                <ul class="list-unstyled mb-0">
                  <li v-for="item in sub.items" :key="item.name" class="mb-3">
                    <a :href="item.url" target="_blank" rel="noopener noreferrer" class="text-decoration-none fw-medium"
                      style="color: var(--app-primary);">
                      {{ item.name }}
                    </a>
                    <p style="color: var(--app-text-secondary); margin-bottom: 0; margin-left: 0.75rem;">
                      {{ item.desc }}
                    </p>
                  </li>
                </ul>
              </div>
            </div>
          </section>
        </div>

        <div class="text-center mt-5 pt-4" style="border-top: 1px solid var(--app-custom-border-color);">
          <p style="color: var(--app-text-secondary); margin-bottom: 0;">
            <i class="bi bi-info-circle me-1"></i>
            {{ footerText }}
          </p>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import { useI18n } from 'vue-i18n'
import { loadResources } from '@/utils/contentLoader'

/*
  ResourceView
  - 资源页面
*/
export default {
  name: 'ResourceView',
  setup() {
    const { t, locale } = useI18n()
    return { t, locale }
  },
  data() {
    return { resources: [] }
  },
  computed: {
    pageTitle() {
      return this.t('resources')
    },
    pageSubtitle() {
      return this.t('resourceSubtitle')
    },
    footerText() {
      return this.t('updating')
    }
  },
  created() {
    this.loadResourcesData()
  },
  watch: {
    locale() {
      this.loadResourcesData()
    }
  },
  methods: {
    async loadResourcesData() {
      try {
        this.resources = await loadResources() || [];
      } catch (error) {
        console.error('Failed to load resources data:', error);
        this.resources = [];
      }
    }
  }
}
</script>

<style scoped>
.typography-body {
  font-size: 1.125rem;
  line-height: 1.8;
}

.typography-body p {
  margin-bottom: 0.75rem;
}
</style>
